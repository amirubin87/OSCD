package WithLevels;

public class Edge {
	int lowNodeId;
	int highNodeId;
	double weight;
	
	public Edge(int v, int u, double w){
		lowNodeId = Math.min(v, u);
		highNodeId = Math.max(v, u);
		weight = w;
	}

}
