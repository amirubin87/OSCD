
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

public class SCDGraphMetaData {
	UndirectedUnweightedGraph g;    
    public Map<Integer, Long> T;
    public Map<Integer, Set<Integer>> VT;  
    public Map<Integer, Map<Integer, Integer>> Intersection_c1_c2;        
    public Map<Integer, Set<Integer>> com2nodes;
    public Map<Integer, Set<Integer>> node2coms;
    
    //cons
    public SCDGraphMetaData(){
    	Intersection_c1_c2 = new HashMap<Integer, Map<Integer, Integer>>();
    	T = new HashMap<Integer, Long>();
    	VT = new HashMap<Integer, Set<Integer>>();
    	com2nodes = new HashMap<Integer, Set<Integer>>();
    	node2coms = new HashMap<Integer, Set<Integer>>();
    }
    
    public SCDGraphMetaData(UndirectedUnweightedGraph graph){        
        //VerifyNodesNumbers(graph);
    	this(); 
    	g = graph;
    	T = graph.Triangles();
    	VT = graph.VTriangles();
    	Integer count = 0;
    
        for (Integer node : graph.nodes()){
            Intersection_c1_c2.put(count, new HashMap<Integer, Integer>());
            Set<Integer> comm = new HashSet<Integer>();
            comm.add(node);
            com2nodes.put(count, comm);
            Set<Integer> commId = new HashSet<Integer>();
            commId.add(count);
            node2coms.put(node, commId);          
            count++;
        }

    }

    public SCDGraphMetaData(UndirectedUnweightedGraph graph, Map<Integer,Set<Integer>> comms){
    	this(); 
    	g = graph;
    	T = graph.Triangles();   
    	VT = graph.VTriangles(); 	
        for (Entry<Integer, Set<Integer>> comm :comms.entrySet()){
        	int commID =comm.getKey();
        	Set<Integer> nodes =comm.getValue();        
            com2nodes.put(commID,nodes);
            Intersection_c1_c2.put(commID,new HashMap<>());
            for (int node : nodes){
            		Set<Integer> commInSet = new HashSet<>();
            		commInSet.add(commID);
	                node2coms.put(node,commInSet);
            }
        }
    }
    
    public SCDGraphMetaData(UndirectedUnweightedGraph graph, Map<Integer,Set<Integer>> comms, boolean partitionIsFromFile){
    	this(); 
    	g = graph;
    	T = graph.Triangles();   
    	VT = graph.VTriangles(); 
    	int maxCommIDSeen = 0;
    	int commID = 0;
        for (Entry<Integer, Set<Integer>> comm :comms.entrySet()){
        	commID =comm.getKey();
        	maxCommIDSeen = Math.max(maxCommIDSeen, commID);
        	Set<Integer> nodes =comm.getValue();        
            com2nodes.put(commID,nodes);
            Intersection_c1_c2.put(commID,new HashMap<>());
            for (int node : nodes){
            		Set<Integer> commsNodeIsIn = node2coms.get(node);
            		if(commsNodeIsIn == null){
	            		commsNodeIsIn = new HashSet<>();	            		
		                node2coms.put(node,commsNodeIsIn);
            		}
        			commsNodeIsIn.add(commID);            		
            }
        }
        
        if(partitionIsFromFile){
        	for( Entry<Integer, Set<Integer>> AcommIdAndNodes: com2nodes.entrySet()){
        		for( Entry<Integer, Set<Integer>> BcommIdAndNodes: com2nodes.entrySet()){
        			int AcommId = AcommIdAndNodes.getKey();
        			int BcommId = BcommIdAndNodes.getKey();
        			if (AcommId < BcommId){
        				int intersectionSize = Utills.IntersectionSize(AcommIdAndNodes.getValue(), BcommIdAndNodes.getValue());
        				Map<Integer, Integer> AcommDictionary = Intersection_c1_c2.get(AcommId);
        				AcommDictionary.put(BcommId, intersectionSize);
        			}
        		}
        	}
        	commID = maxCommIDSeen +1;
        	for(int node : g.nodes()){
        		Set<Integer> nodesComms = node2coms.get(node);
        		if(nodesComms == null){
        			nodesComms = new HashSet<>();
        			nodesComms.add(commID);
        			node2coms.put(node,nodesComms);
        			Set<Integer> nodeInSet = new HashSet<>();
        			nodeInSet.add(node);
        			com2nodes.put(commID, nodeInSet);
        			commID++;
        		}
        	}
        }
    }
    
	public SCDGraphMetaData(SCDGraphMetaData ORIGINALmetaData) {
    	g=ORIGINALmetaData.g;
    	T=Utills.CopyMapIntLong(ORIGINALmetaData.T);
    	VT=Utills.CopyMapIntSet(ORIGINALmetaData.VT);
    	Intersection_c1_c2 = Utills.CopyMapIntMapIntInt(ORIGINALmetaData.Intersection_c1_c2);
    	com2nodes = Utills.CopyMapIntSet(ORIGINALmetaData.com2nodes);
    	node2coms = Utills.CopyMapIntSet(ORIGINALmetaData.node2coms);
	}

	public void ClearCommsOfNode(Integer node){
    	Set<Integer> commsSet = node2coms.get(node); 
    	
    	//update intersection ratio
    	UpdateIntersectionRatioRemove(commsSet);   
    	
    	//Symbolic change
        for(Integer comm : commsSet){
            com2nodes.get(comm).remove(node);
        }

        node2coms.put(node, new HashSet<Integer>());
    }
   
	private void UpdateIntersectionRatioRemove(Set<Integer> c_v) {
		// Intersection size update
        Integer[] comms = new Integer[c_v.size()];
        int k = 0;
        for(Integer comm : c_v){
        	comms[k] = comm;
        	k++;      	
        }           
        
        for (int i = 0; i <comms.length ; i++){
        	for (int j = i+1; j < comms.length ; j++){
        		int x = comms[i];
        		int y = comms[j];
        		
        		Integer min = Math.min(x, y);
        		Integer max = Math.max(x, y);
        		Intersection_c1_c2.get(min).put(max, Intersection_c1_c2.get(min).get(max)-1);
        	}
        }		
	}
	
	private void UpdateIntersectionRatioAdd(Set<Integer> c_v) {
		// Update Intersection_c1_c2
        Integer[] commsArray = new Integer[c_v.size()];
        int k = 0;
        for(Integer comm : c_v){
        	commsArray[k] = comm;
        	k++;      	
        }
        
	    for (int i = 0; i <commsArray.length ; i++){
	    	for (int j = i+1; j < commsArray.length ; j++){
	    		int x = commsArray[i];
	    		int y = commsArray[j];
	    		Integer lowComm = Math.min(x, y);
	    		Integer highComm = Math.max(x, y);	
		        Map<Integer,Integer> lowCommDic = Intersection_c1_c2.get(lowComm);
			    if ( lowCommDic == null){
			    	lowCommDic = new HashMap<Integer, Integer>();
		            Intersection_c1_c2.put(lowComm, lowCommDic);
			    }	    
			    
			    Integer intersectionSize = lowCommDic.get(highComm);
		        if ( intersectionSize == null){
		        	intersectionSize = 0;
		        }   
		        
		        intersectionSize++;
		    
		        lowCommDic.put(highComm, intersectionSize);
	    	}
	    }
	}

	public Map<Integer[],Double> SetCommsForNode(Integer node, Set<Integer> comms, boolean shouldMergeComms){		
    	
		UpdateIntersectionRatioAdd(comms);
		
		Map<Integer[],Double> commsCouplesIntersectionRatio = new HashMap<Integer[],Double>();
	    
		// Find intersection ration for merge
        Integer[] commsArray = new Integer[comms.size()];
        int k = 0;
        for(Integer comm : comms){
        	commsArray[k] = comm;
        	k++;      	
        }        
        
		// Symbolic add
	    for (Integer comm : comms){
	        com2nodes.get(comm).add(node);	        
	    }
	    
	    node2coms.put(node, comms);
        if(shouldMergeComms){
		    for (int i = 0; i <commsArray.length ; i ++){
		    	for (int j = i+1; j < commsArray.length ; j++){
		    		int x = commsArray[i];
		    		int y = commsArray[j];
		    		Integer lowComm = Math.min(x, y);
		    		Integer highComm = Math.max(x, y);
			        double intersectionRatio = (double)Intersection_c1_c2.get(lowComm).get(highComm)/(double)Math.min(com2nodes.get(lowComm).size(), com2nodes.get(highComm).size());
			        Integer[] sortedComms= new Integer[]{lowComm,highComm};
			        commsCouplesIntersectionRatio.put(sortedComms, intersectionRatio);
		    	}
		    }
        }
	    

	    return commsCouplesIntersectionRatio;
    }

	public void AddCommForNode(int node, int c1){
        
		Set<Integer> comms = node2coms.get(node);
        if(comms.contains(c1)){
        	return;
        }
        
        for (int comm : comms){
            int lowComm = Math.min(comm, c1);
            int highComm = Math.max(comm, c1);

            Map<Integer,Integer> lowCommDic = Intersection_c1_c2.get(lowComm);
		    if ( lowCommDic == null){
		    	lowCommDic = new HashMap<Integer, Integer>();
	            Intersection_c1_c2.put(lowComm, lowCommDic);
		    }	    
		    
		    Integer intersectionSize = lowCommDic.get(highComm);
	        if ( intersectionSize == null){
	        	intersectionSize = 0;
	        }   
	        
	        intersectionSize++;
	        lowCommDic.put(highComm,intersectionSize);
        }
        node2coms.get(node).add(c1);
        com2nodes.get(c1).add(node);
	}
	
	public void RemoveCommForNode(int node, int c1){        
        if(!node2coms.get(node).contains(c1)){
        	return;
        }        
        
        node2coms.get(node).remove(c1);
        com2nodes.get(c1).remove(node);
        Set<Integer> comms = node2coms.get(node);
		
        for (int comm : comms){
            int lowComm = Math.min(comm, c1);
            int highComm = Math.max(comm, c1);

            Map<Integer,Integer> lowCommDic = Intersection_c1_c2.get(lowComm);
		    if ( lowCommDic != null){
		    	Integer intersectionSize = lowCommDic.get(highComm);
		    	if ( intersectionSize != null){
		    		intersectionSize--;
		    		lowCommDic.put(highComm,intersectionSize);
		    		}
		    }
        }    
        
	}
                    
	public void SymbolicClearComm(Integer comm){    	
    	Set<Integer> nodes = new HashSet<Integer>(com2nodes.get(comm));
    	for (Integer node: nodes){
    		SymbolicRemoveNodeFromComm(node, comm);
    	}
    }
    
	public void SymbolicRemoveNodeFromComm(Integer node,Integer comm){
        node2coms.get(node).remove(comm);
        com2nodes.get(comm).remove(node);        
    }
    
	public void SymbolicAddNodeToComm(Integer node, Integer comm) {
		node2coms.get(node).add(comm);
		com2nodes.get(comm).add(node);		
	}

}
    
/*
 * 
__a
uthor__ = 't-amirub'
import networkx as nx
import itertools
import copy

class SCDMetaData :
    T = dict()
    com2nodes= dict()
    node2coms = dict()
    Intersection_c1_c2 = dict()


    def __str__(self) :
        return ( " T : " + str(self.T)
            + " com2nodes : " + str(self.com2nodes)
            + " node2coms : " + str(self.node2coms)
            + " Intersection_c1_c2 : " + str(self.Intersection_c1_c2))


   def copy(self) :
        new_MetaData = SCDMetaData()
        new_MetaData.T = copy.deepcopy(self.T)
        new_MetaData.com2nodes = copy.deepcopy(self.com2nodes)
        new_MetaData.node2coms = copy.deepcopy(self.node2coms)
        new_MetaData.Intersection_c1_c2 = self.Intersection_c1_c2.copy()
        return new_MetaData

    def init(self, graph, initPart):
        """Initialize the SCDMetaData of a graph with every node in one community"""
        self.T = nx.triangles(graph)
        self.com2nodes = dict()
        self.node2coms = dict()
        self.Intersection_c1_c2 = dict()
        for comm, nodes in initPart.items():
            self.com2nodes[comm] = nodes
            self.Intersection_c1_c2[comm] = dict([])
            for node in nodes:
                self.node2coms[node] = set([comm])


    def ClearComms(self, node):
        comms = self.node2coms.get(node)
        if comms ==None:
            comms = set()
        for comm in comms:
            self.com2nodes[comm].remove(node)
        self.node2coms[node] = set()
        for pairComms in itertools.combinations(comms,2):
            self.Intersection_c1_c2[min(pairComms)][max(pairComms)] -= 1


    def AddCommForNode(self, node, c1):
        comms = self.node2coms[node]
        for comm in comms:
            lowComm = min(comm, c1)
            highComm = max(comm, c1)

            if self.Intersection_c1_c2.get(lowComm) == None:
                self.Intersection_c1_c2[lowComm] = dict([])
            if self.Intersection_c1_c2[lowComm].get(highComm) == None:
                self.Intersection_c1_c2[lowComm][highComm] = 0
            self.Intersection_c1_c2[lowComm][highComm] += 1
        self.node2coms[node].add(c1)
        self.com2nodes[c1].add(node)


    def RemoveCommForNode(self, node,c1):
        if self.node2coms.get(node) != None:
            self.node2coms[node].discard(c1)
        if self.com2nodes.get(c1) != None:

            self.com2nodes[c1].remove(node)
        for comm in self.node2coms[node]:
            if self.Intersection_c1_c2.get(min(comm, c1)) != None:
                if self.Intersection_c1_c2[min(comm, c1)].get(max(comm,c1)) !=None:
                    self.Intersection_c1_c2[min(comm, c1)][max(comm,c1)] -= 1

    def SetCommsForNode(self, node, comms):
        commsCouplesIntersectionRatio = []

        for pairComms in itertools.combinations(comms, 2):
            lowComm = min(pairComms)
            highComm = max(pairComms)
            if self.Intersection_c1_c2.get(lowComm) == None:
                self.Intersection_c1_c2[lowComm] = dict([])
            if self.Intersection_c1_c2[lowComm].get(highComm) == None:
                self.Intersection_c1_c2[lowComm][highComm] = 0
            intersectionSize = self.Intersection_c1_c2[lowComm][highComm] + 1
            #self.Intersection_c1_c2[lowComm][highComm] = intersectionSize will be increased in AddCommForNode
            intersectionRatio = intersectionSize/min(len(self.com2nodes[lowComm]), len(self.com2nodes[highComm]))
            commsCouplesIntersectionRatio.append((lowComm, highComm, intersectionRatio))
        for comm in comms:
            self.AddCommForNode(node, comm)
        return commsCouplesIntersectionRatio

                                                                                                         93,0-1        Bot
*/