package WithLevels;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;


public class UndirectedGraph {
	private Set<Integer> nodes;
	private Map<Integer,Set<Integer>> neighbors;
	private double weight_of_edges;
	private double num_of_edges;
	private Map<Integer, Double> T;
	private Map<Integer, Map<Integer, Double>> EdgesWeights;
	private Map<Integer, Set<Integer>> VT;
	private int MaxNodeId=0;
	private Map<Integer, Double> degrees;
	private Map<Integer, Double> VTWeights;
	
	// CONS
	public UndirectedGraph(){
		nodes = new HashSet<>();
		neighbors = new HashMap<>();
		weight_of_edges = 0.0;
		num_of_edges = 0;
		T = new HashMap<>();
		EdgesWeights = new HashMap<>();
		VT = new HashMap<>();
		MaxNodeId=0;
		degrees = new HashMap<>();
		VTWeights = new HashMap<>();
	}
	
	public UndirectedGraph(Path p) throws IOException{
		degrees = new HashMap<>();
		EdgesWeights = new HashMap<>();
		nodes = new HashSet<Integer>();
		neighbors = new HashMap<Integer,Set<Integer>>();
		List<String> lines= Files.readAllLines(p);
		int v,u,max,min;
		double weight;
		Set<Integer> vNeig, uNeig;
		for (String line : lines){
			String[] parts = line.split(" |\t");
			v = Integer.parseInt(parts[0].trim());
			u = Integer.parseInt(parts[1].trim());
			if(parts.length>2){
				weight = Double.parseDouble(parts[2].trim());
			}
			else{
				weight=1.0;
			}
			
			max = Math.max(u, v);
			min = Math.min(u, v);
			if(max>MaxNodeId){
				MaxNodeId=max;
			}
			
			vNeig= neighbors.get(v);
			if(vNeig == null){
				vNeig = new HashSet<Integer>();
				neighbors.put(v, vNeig);
			}
			uNeig= neighbors.get(u);
			if(uNeig == null){
				uNeig = new HashSet<Integer>();
				neighbors.put(u, uNeig);
			}
			if(v!=u && !vNeig.contains(u)){
				vNeig.add(u);
				uNeig.add(v);
				weight_of_edges= weight_of_edges+weight;
				num_of_edges++;
				nodes.add(v);
				nodes.add(u);	
				AddEdge(min,max,weight);
			}			
		}
		CalcTrianglesAndVT();
	}
	
	// Getters
	
	public int GetNumberOfNodes() {		
		return nodes.size();
	}

	public Set<Integer> GetNodes() {		
		return nodes;
	}
		 
	public double GetWeightOfAllEdges() {		
		return weight_of_edges;
	}
	
	public Map<Integer, Double> GetTriangles() {
		return T;
	}
	
	public Map<Integer, Set<Integer>> GetVTriangles() {
		return VT;
	}
	
	public int GetMaxNodeId() {		
		return MaxNodeId;
	}
	
	//Private
	private void AddEdge(int smallNode, int bigNode, double weight) {
		Map<Integer,Double> smallNodeMap = EdgesWeights.get(smallNode);
		if(smallNodeMap == null){
			smallNodeMap = new HashMap<>();
			EdgesWeights.put(smallNode,smallNodeMap);
		}
		smallNodeMap.put(bigNode, weight);
		
		Double oldDegree = degrees.get(smallNode);
		if(oldDegree == null){
			degrees.put(smallNode, weight);
		}
		else{
			degrees.put(smallNode, weight + oldDegree);
		}
		
		oldDegree = degrees.get(bigNode);
		if(oldDegree == null){
			degrees.put(bigNode, weight);
		}
		else{
			degrees.put(bigNode, weight + oldDegree);
		}
	}
	
	private double ClustringPerNode(Integer node) {
		double d = NodeDegreeWeight(node);
		if (d <1){
			return 0;
		}		
		return (double)(2*T.get(node))/(d*(d-1));
	}
		
   private void CalcTrianglesAndVT() {
    	T = new HashMap<Integer, Double>();    	 
    	VT = new HashMap<Integer, Set<Integer>>();
    	VTWeights = new HashMap<>();
    	for(int v : nodes){    		
    		T.put(v,(double) 0);    		
    		VTWeights.put(v,(double) 0);
    		VT.put(v,new HashSet<>());
    	}
    	
    	Set<Integer> vTriangle, uTriangle, wTriangle;    	
    	for(int v : nodes){
    		Set<Integer> vNeighbors = neighbors(v);
    		vTriangle = VT.get(v);
    		for( int u : vNeighbors){
    			if(u > v){
    				uTriangle = VT.get(u);
    				for(int w : neighbors(u)){
    					if (w > u && vNeighbors.contains(w)){
    						wTriangle = VT.get(w);
    						AddNodeToVTAndUpdateVTWeigh(v, vTriangle,u);
    						AddNodeToVTAndUpdateVTWeigh(v, vTriangle,w);
    						AddNodeToVTAndUpdateVTWeigh(u, uTriangle,v);
    						AddNodeToVTAndUpdateVTWeigh(u, uTriangle,w);
    						AddNodeToVTAndUpdateVTWeigh(w, wTriangle,v);
    						AddNodeToVTAndUpdateVTWeigh(w, wTriangle,u);
    						double triangleW = TriangleWeight(v,u,w);
    						T.put(v, T.get(v)+triangleW);
    						T.put(u, T.get(u)+triangleW);
    						T.put(w, T.get(w)+triangleW);
    						
    					}
    				}
    			}
    		}
    	}
	}

   private void AddNodeToVTAndUpdateVTWeigh(int v, Set<Integer> vTriangle, int u) {
	   int sizeBeforeAdd = vTriangle.size();
	   vTriangle.add(u);
	   if(sizeBeforeAdd!=vTriangle.size()){
		   VTWeights.put(v, VTWeights.get(v) + EdgeWeight(u, v));
	   }
	
}

   private double TriangleWeight(int v, int u, int w) {
		return Math.min(EdgeWeight(v, u), Math.min(EdgeWeight(v, w), EdgeWeight(u, w)));
	}

	
	//public
   
	public Map<Integer,Double> Clustring() {
		Map<Integer,Double> ans = new HashMap<>();
		for(int node:nodes){
			ans.put(node, ClustringPerNode(node));
		}
		return ans;
	}

	public double EdgeWeight(int u, int v) {
		int smallNode=Math.min(u, v);
		int bigNode = Math.max(u, v);
		Map<Integer,Double> smallNodeMap = EdgesWeights.get(smallNode);
		if(smallNodeMap == null){
			return 0;
		}
		Double ans =  smallNodeMap.get(bigNode);
		if(ans == null){
			return 0;
		}
		return ans;
	}

	public String toString() {		
		return "Num of nodes: " + nodes.size() + " . Weight of edges: " + weight_of_edges;
	}
	
	public Set<Integer> neighbors(int node) {		
		return neighbors.get(node);
	}
	
	public double NodeDegreeWeight(Integer node) {		
		return degrees.get(node);
	}

	
	public Map<Integer, Double> GetVTrianglesWeights() {		
		return VTWeights;
	}

	public double GetAverageWeight() {
		return weight_of_edges/num_of_edges;
	}

	private void AddNode(int nodeId) {
		if(!nodes.contains(nodeId)){
			nodes.add(nodeId);
			if(MaxNodeId < nodeId ){
				MaxNodeId = nodeId;
			}
			neighbors.put(nodeId, new HashSet<>());
			EdgesWeights.put(nodeId, new HashMap<>());
			degrees.put(nodeId, (double) 0);
		}
	}

	public void AddEdges(int commId, Map<Integer, Double> edgesToadd) {
		AddNode(commId);
		for(Entry<Integer, Double> neighWeight : edgesToadd.entrySet()){
			Integer neighbor= neighWeight.getKey();
			AddNode(neighbor);
			double weight= neighWeight.getValue();
			weight_of_edges= weight_of_edges+weight;
			num_of_edges++;
			neighbors.get(commId).add(neighbor);
			neighbors.get(neighbor).add(commId);
			int min = Math.min(neighbor,commId);
			int max = Math.max(neighbor,commId);
			AddEdge(min, max, weight);
			
		}
			
	}

	public void CalcGraphMetrics() {
		CalcTrianglesAndVT();		
	}

	
}

