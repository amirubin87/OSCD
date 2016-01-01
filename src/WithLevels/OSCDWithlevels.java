package WithLevels;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.text.html.parser.Entity;

import java.util.Map.Entry;

public class OSCDWithlevels {
	/***************************    Fields  ************************************/
	// DS
	public List<UndirectedGraph> LevelsGraphs;
	public UndirectedGraph OriginalG;
	public UndirectedGraph g;
	
	public List<SCDWithLevelsGraphMetaData> LevelsMetaData;
	public SCDWithLevelsGraphMetaData OriginalMetaData;
	public SCDWithLevelsGraphMetaData metaData;
	
	// params
	public double[] betas;
	public double alpha;
	public String outputPath;
	public int iteratioNumToStartMerge;
	public int maxIterationsToRun;
	public int percentageOfStableNodes;
	public String pathToGraph;	
	
	// results
	public List<Map<Integer,Set<Integer>>> Dendogram;
	
	/***************************    CONS ************************************/
	 
	public OSCDWithlevels(String pathToGraph, double[]betas, double alpha, String outputPath, int iteratioNumToStartMerge, int maxIterationsToRun, int percentageOfStableNodes) throws IOException{
		this.percentageOfStableNodes= percentageOfStableNodes;
		this.betas= betas;
		this.alpha = alpha;
		this.outputPath =outputPath;
		this.iteratioNumToStartMerge = iteratioNumToStartMerge;
		this.maxIterationsToRun = maxIterationsToRun;
		this.pathToGraph = pathToGraph;
		
		this.OriginalG = new UndirectedGraph(Paths.get(pathToGraph));		
		
		Map<Integer, Set<Integer>> firstPart = GetFirstPartition(OriginalG);		
		this.OriginalMetaData = new SCDWithLevelsGraphMetaData(OriginalG,firstPart,false);
	}
	
	public OSCDWithlevels(String pathToGraph, String pathToPartition, double[]betas, double alpha, String outputPath, int iteratioNumToStartMerge, int maxIterationsToRun, int percentageOfStableNodes) throws IOException{
		this.percentageOfStableNodes= percentageOfStableNodes;
		this.betas= betas;
		this.alpha = alpha;
		this.outputPath =outputPath;
		this.iteratioNumToStartMerge = iteratioNumToStartMerge;
		this.maxIterationsToRun = maxIterationsToRun;
		this.pathToGraph = pathToGraph;
		
		this.OriginalG = new UndirectedGraph(Paths.get(pathToGraph));
		
		Map<Integer, Set<Integer>> firstPart = GetPartitionFromFile(pathToPartition);		
		this.OriginalMetaData = new SCDWithLevelsGraphMetaData(OriginalG,firstPart, false);
	}
	
	
	
	/***************************    find comms  ************************************/
	
	public void FindCommunities() throws IOException{
		this.Dendogram = new ArrayList<>();
		for (double betta : betas){
			System.out.println("");
			System.out.println("                       Input: " + pathToGraph);
			System.out.println("                       betta: " + betta);
			// Create a copy of the original meta data
			this.LevelsMetaData = new ArrayList<>();		
			LevelsMetaData.add(new SCDWithLevelsGraphMetaData(OriginalMetaData));
			this.LevelsGraphs = new ArrayList<>();
			LevelsGraphs.add(OriginalG);
			int level = 0;
			// We break when no comms are found, like done in OSLOM
			while(true){
				metaData = new SCDWithLevelsGraphMetaData(LevelsMetaData.get(level));
				g= LevelsGraphs.get(level);
				Map<Integer,Set<Integer>> levelComms = FindCommunities(betta,level);
				//we stop once we have less than 100 communities
				int nonEmpty = CountNonEmpty(levelComms);
				if(nonEmpty <100){
					break;
				}
				//Store in dendogram
				AddToDendogram(level,levelComms);
				// output level				
				WriteToFile(Dendogram.get(level), betta, level);
				if(nonEmpty <2){
					break;
				}
				//Create next level graph and meta data
				CreateNextLevel();
				level++;
			}
		}
	}	
	

	private int CountNonEmpty(Map<Integer, Set<Integer>> levelComms) {
		int ans= 0;
		for(  Set<Integer> comm : levelComms.values()){
			if(comm.size()>0)
				ans++;
		}
		return ans;
	}

	// In every level all nodes will appear!
	public void AddToDendogram(int level, Map<Integer, Set<Integer>> levelComms) {
		if (level ==0){
			Dendogram.add(level, levelComms);
			return;
		}
		Map<Integer, Set<Integer>> currentLevel = new HashMap<>();
		Map<Integer, Set<Integer>> lowerLevel = Dendogram.get(level-1);
		
		for( Entry<Integer, Set<Integer>> topLevelComm: levelComms.entrySet()){
			int commID = topLevelComm.getKey();

			Set<Integer> nodes = new HashSet<>();
			Set<Integer> topLevelNodes = topLevelComm.getValue();

			for(int node:topLevelNodes ){

				nodes.addAll(lowerLevel.get(node));
			}
			currentLevel.put(commID, nodes);
		}		
		Dendogram.add(level, currentLevel);
	}

	private Map<Integer,Set<Integer>> FindCommunities(double betta, int level) {
	    int numOfStableNodes = 0;
	    int amountOfScans = 0;
	    int n = g.GetNumberOfNodes();
	    int tenPercent = n/10+1;
	    int numOfStableNodesToReach = n*percentageOfStableNodes/100;
	    System.out.println();
	    System.out.println();
	    while (numOfStableNodes < numOfStableNodesToReach && amountOfScans < maxIterationsToRun){	    	
	    	System.out.println();
	    	System.out.println("Input: " +pathToGraph + " betta: " + betta + "            Num of iter: " + amountOfScans);
	    	System.out.println("Amount Of Stable Nodes: " + numOfStableNodes);
	    	System.out.print("progress in current iteration: ");
	    	numOfStableNodes=0;
	        amountOfScans++;
	        int nodeCounter= 0;
	        for (Integer node : g.GetNodes()){
	        	nodeCounter++;	        	
	        	if ((nodeCounter%tenPercent) == 0){
	        		System.out.print(nodeCounter/tenPercent*10 + "%  ");
	        	}
	            Set<Integer> c_v_original = metaData.node2coms.get(node);	            
	            metaData.ClearCommsOfNode(node);
	            Map<Integer, Double> comms_inc = new HashMap<Integer, Double>();
	            Set<Integer> neighborComms = Find_Neighbor_Comms(node);
	            for (Integer neighborComm : neighborComms){
	                double inc= Calc_WCC(neighborComm, node);
	                comms_inc.put(neighborComm, inc);
	            }	            
	            Set<Integer> c_v_new =Keep_Best_Communities(comms_inc, betta);
	           
	            boolean shouldMergeComms = amountOfScans>iteratioNumToStartMerge;
				Map<Integer[],Double> commsCouplesIntersectionRatio = metaData.SetCommsForNode(node, c_v_new, shouldMergeComms );
	            boolean haveMergedComms = false;
	            if(shouldMergeComms){
	            	haveMergedComms = FindAndMergeComms(commsCouplesIntersectionRatio);
	            }            
	            
	            if (!haveMergedComms && c_v_new.equals(c_v_original)){
	            	numOfStableNodes++;
	            }
	        }
        }    
	    if (amountOfScans >= maxIterationsToRun){
	        System.out.println(String.format("NOTICE - THE ALGORITHM HASNT STABLED. IT STOPPED AFTER SCANNING ALL NODES FOR %1$d TIMES.",maxIterationsToRun));
	    }
	    return metaData.com2nodes;
	}	  
		
	private boolean FindAndMergeComms (Map<Integer[],Double> commsCouplesIntersectionRatio){
	    boolean haveMergedComms = false;
	    //Set<Integer> commsToClean = new HashSet<Integer>();
	    for (Entry<Integer[],Double > c1c2intersectionRate : commsCouplesIntersectionRatio.entrySet()){	    	
	        if(c1c2intersectionRate.getValue()>alpha){
	        	Integer[] c1c2 = c1c2intersectionRate.getKey();
	        	//commsToClean.add(c1c2[0]);
	        	MergeComms(c1c2);
	        	haveMergedComms = true;
	        }
	    }
	    //ClearNodesFromComms(commsToClean);
	    return haveMergedComms;
	}

	private void MergeComms(Integer[] commsToMerge){
		Integer c1 = commsToMerge[0];
		Integer c2 = commsToMerge[1];
		List<Integer> copyOfC1= new ArrayList<>(metaData.com2nodes.get(c1));
		List<Integer> copyOfC2= new ArrayList<>(metaData.com2nodes.get(c2));
	    for (Integer node : copyOfC1){	 
	    	metaData.RemoveCommForNode(node,c1);
	        if(!copyOfC2.contains(node)){	        	
	        	metaData.AddCommForNode(node,c2);
	        }	        
	    }
	}
	
	private Set<Integer> Keep_Best_Communities(Map<Integer, Double>comms_imps,double betta){
	    double bestImp = 0;
	    for( double imp : comms_imps.values()){
	    	bestImp = Math.max(bestImp, imp);
	    }
	    Set<Integer> bestComs = new HashSet<Integer>();
	    for(Entry<Integer, Double> entry: comms_imps.entrySet()){
	    		 if (entry.getValue()*betta > bestImp){
	    				 bestComs.add(entry.getKey());
	    		 }
	    }
	    return bestComs;
	}	

	private Set<Integer> Find_Neighbor_Comms(Integer node){
	    Set<Integer>neighborComms = new HashSet<Integer>();
	    for (Integer neighbor : g.neighbors(node)){
	        neighborComms.addAll(metaData.node2coms.get(neighbor));
	    }
    return neighborComms;
    }
	
	private void WriteToFile(Map<Integer, Set<Integer>> comms, double betta, int level) throws FileNotFoundException, UnsupportedEncodingException {
		PrintWriter writer = new PrintWriter(outputPath + betta + "-L-"+ level + ".txt", "UTF-8");
		for ( Set<Integer> listOfNodes : comms.values()){
			if(listOfNodes.size()>2){
				for(int node : listOfNodes){
					writer.print(node + " ");
				}
				writer.println("");
			}
		}		
		writer.close();	
	}

	public static Map<Integer,Set<Integer>> GetFirstPartition(UndirectedGraph G){
		Map<Integer,Set<Integer>> result = new HashMap<>();
		Map<Integer, Double> CC = G.Clustring();		
	    Map<Integer, Double> sorted_CC = MapUtil.sortByValue(CC);
	    double maxSeenSoFar=1.0;    
	    boolean[] isVisited = new boolean[G.GetMaxNodeId()+1];	    
	    int commID=0;
	    int numOfNodes = sorted_CC.size();
		int tenPercent = numOfNodes/10+1;
		int nodeCounter = 0;
	    System.out.println();
		System.out.print("Finding first partitioning. Progress: ");
	    for (int v : sorted_CC.keySet()){
	    	nodeCounter++;
			if ((nodeCounter%tenPercent) == 0){
        		System.out.print(nodeCounter/tenPercent*10 + "%  ");
        	}
	    	if(maxSeenSoFar<CC.get(v)){
	    		throw(new RuntimeException(String.format("sortedCC was not sorted. node: %1$d.", v)));
	    	}	    	
	        if (!isVisited[v]){
	            isVisited[v]= true;
	            Set<Integer> vSet = new HashSet<>();
	            vSet.add(v);
	            result.put(commID, vSet);
	            for(int  neigh : G.neighbors(v)){
	            	if (!isVisited[neigh]){
	            		isVisited[neigh]= true;
	                    result.get(commID).add(neigh);
	                }
	            }
	            commID+=1;
	        }
	    }
	    
	    return result;
	}
	
	public static Map<Integer,Set<Integer>> GetFirstPartitionCliques(UndirectedGraph G){
		Set<Integer> hasComm = new HashSet<>();
		boolean vHasComm=false;
		Map<Integer,Set<Integer>> result = new HashMap<>();
			    
	    int commID=0;
	    Set<Integer> nodes = G.GetNodes();
	    for (int v :nodes){
	    	if(hasComm.contains(v)){
	    		continue;
	    	}
	    	vHasComm = false;
	    	Set<Integer> vNeigh = G.neighbors(v);
	    	for(int u:vNeigh){
	    		if(vHasComm){
	    			break;
	    		}
	    		if(!hasComm.contains(u)){
	    			Set<Integer> UVNeigh= Utills.Intersection(vNeigh, G.neighbors(u));
	    			for(int w:UVNeigh){
	    				if(vHasComm){
	    	    			break;
	    	    		}
	    				if(!hasComm.contains(w)){	
	    					for(int z : Utills.Intersection(UVNeigh, G.neighbors(w))){
	    						if(vHasComm){
	    			    			break;
	    			    		}
	    						if(!hasComm.contains(z)){	    					
			    					Set<Integer> comm = new HashSet<>();
			    					comm.add(v);
			    					comm.add(u);
			    					comm.add(w);
			    					comm.add(z);
			    					result.put(commID, comm);
			    					commID+=1;
			    					hasComm.add(v);
			    					hasComm.add(u);
			    					hasComm.add(w);
			    					hasComm.add(z);
			    					vHasComm = true;
			    					break;
			    				}	    						
	    					}
	    				}
	    			}
	    		}
	    	}
	    	if(!vHasComm){
	    		Set<Integer> comm = new HashSet<>();
				comm.add(v);
	    		result.put(commID, comm);
	    		commID+=1;
	    		hasComm.add(v);
	    	}
	    }
	    return result;
	}
	
	public double OLD_WCC(int comm, int  node){
		int n =(int) g.GetNumberOfNodes(); 
		Set<Integer> commMembers = metaData.com2nodes.get(comm);
		 double t = calcT(commMembers, node);
		 double divesor = commMembers.size() + n*((metaData.T.get(node) - t));
		 if (divesor==0)
		        return 0;
		 return (double)n*t/(double)divesor;
	}
				   
	public double Calc_WCC(int comm, int  x){	    
		Set<Integer> commMembers = metaData.com2nodes.get(comm);
		double TxV = metaData.T.get(x);	    
	    if (TxV==0){
	        return 0;
	    }
	    
		double TxC = calcT(commMembers, x);	    
		if(TxC == 0){
			return 0;
		}
		BigDecimal partA = new BigDecimal(TxC).divide(new BigDecimal(TxV),10, BigDecimal.ROUND_DOWN); 
	    
	    double VTxV = metaData.VTWeight.get(x);
	    if(VTxV == 0){
			return 0;
		}
	    double VTxVnoC = calcVTWithoutComm(commMembers, x);	    
	    double divesor = (double)(commMembers.size()*g.GetAverageWeight() +(VTxVnoC));	    
	    if (divesor==0){
	        return 0;
	    }	    
	    BigDecimal partB = new BigDecimal(VTxV).divide(new BigDecimal(divesor),10, BigDecimal.ROUND_DOWN);	
	    double ans = (partA.multiply(partB)).doubleValue();
	    
	    return ans;
	    
	}

	private double calcVTWithoutComm(Set<Integer> commMembers, int node) {		
		Set<Integer> nodesWithTriangle = metaData.VT.get(node);
		double weightsOfEdgesToNodesWithTriangle = metaData.VTWeight.get(node);
		double weightsOfEdgesToCommMembersWithTriangles = SumWeightsFromNodeToNodes(node, Utills.Intersection(nodesWithTriangle, commMembers));		
		return weightsOfEdgesToNodesWithTriangle - weightsOfEdgesToCommMembersWithTriangles;
	}

	private double SumWeightsFromNodeToNodes(int node, Set<Integer> nodes) {
		double ans = 0;
		for (int v : nodes){
			ans = ans + g.EdgeWeight(node, v);
		}
		return ans;
	}

	private double calcT(Set<Integer> commMembers, int node) {
		double t=0;
	    Set<Integer> neighbours = g.neighbors(node);
	    Set<Integer> neighInComm = Utills.Intersection(commMembers, neighbours);
	    for (int v : neighInComm){
	    	double weightV =g.EdgeWeight(node,v);
	        for (int u : neighInComm){	
	        	double weightUV = g.EdgeWeight(u,v);
	            if (u > v){
	                t = t + Math.min(g.EdgeWeight(u,node), Math.min(weightV,weightUV ));
	                }
	        }
	    }
	    return t;
	}
	
	public Map<Integer,Set<Integer>> GetPartitionFromFile(String partFile) throws IOException{		
		Map<Integer,Set<Integer>> comm2Nodes= new HashMap<Integer,Set<Integer>>();
		List<String> lines= Files.readAllLines(Paths.get(partFile));		
	    int commID=0;
	    for (String line : lines){
	        String[] nodes = line.split(" |\t");	        		 
	    	if (nodes.length >2){
	    		Set<Integer> comm = new HashSet<>();
	    		for (String node : nodes){
	    			comm.add(Integer.parseInt(node.trim()));
	    		}
	            comm2Nodes.put(commID, comm);
	            commID ++;
	    	}
	    }
	    return comm2Nodes;
	}
	
	private void CreateNextLevel(){
		UndirectedGraph nextLevelGraph = new UndirectedGraph();
		for(Entry<Integer, Set<Integer>> comIdNodes: metaData.com2nodes.entrySet()){
			if(comIdNodes.getValue().size()<3){
				continue;
			}
			int commId = comIdNodes.getKey();

			Map<Integer, Double> edgesToadd = FindNeighborCommsAndEdgesWeight(commId, comIdNodes.getValue());

			nextLevelGraph.AddEdges(commId, edgesToadd);
		}
		nextLevelGraph.CalcGraphMetrics();

		LevelsGraphs.add(nextLevelGraph);		
		Map<Integer, Set<Integer>> firstPart = GetFirstPartitionCliques(nextLevelGraph);	
		SCDWithLevelsGraphMetaData nextLevelMetaData = new SCDWithLevelsGraphMetaData(nextLevelGraph, firstPart,true);
		LevelsMetaData.add(nextLevelMetaData);				
	}

	// All will be SMALLER than commID	
	private Map<Integer,Double> FindNeighborCommsAndEdgesWeight(int commId ,Set<Integer> nodesInComm) {
		// go over nodes in comm. find their neigh. Each neighbor- find the comms he is in. If the comm is lower than commId, add to its weight in the map.
		Map<Integer,Double> ans = new HashMap<>();
		Set<Integer> neighbors ;		
		Set<Integer> neighborComms ;
		for (int node :nodesInComm){
			double numOfCommsNodeIsIn = metaData.node2coms.get(node).size(); 
			neighbors = g.neighbors(node);
			for( int neighbor :neighbors){
				neighborComms = metaData.node2coms.get(neighbor);
				double numOfCommsNeighborIsIn = neighborComms.size(); 
				for(Integer neighborComm : neighborComms){
					if(neighborComm < commId && metaData.com2nodes.get(neighborComm).size()>2){
						double w = g.EdgeWeight(node, neighbor)/(numOfCommsNodeIsIn*numOfCommsNeighborIsIn);
						Double weightUntilNow =ans.get(neighborComm);
						if( weightUntilNow == null){
							ans.put(neighborComm, w);
						}
						else{
							ans.put(neighborComm, weightUntilNow + w);
						}
					}
				}
			}
		}
		
		
		return ans;
	}
}





