package WithLevels;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class RunOSCDWithLevels {

	public static void main(String[] args) throws IOException {		
		
		if (false && args.length <3){
			System.out.println("pathToGraph  outputPath  betas  alpha=0.8  iteratioNumToStartMerge=100  maxIterationsToRun=20 percentageOfStableNodes=95");
		}
		else{
			String pathToGraph = "C:/Users/t-amirub/Desktop/amazon/com-amazon.ungraph.txt";
			String outputPath = "C:/Users/t-amirub/Desktop/amazon/SCDLevelsClique/t";
			//String pathToGraph = "C:/Temp/network.txt";
			//String outputPath = "C:/Temp/SCDClique/";
			
			
			double[] betas = {1.1};//1.02,1.01,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,3.5,4.0,4.5}; 
			
			double alpha = 0.8;
			int iteratioNumToStartMerge = 4;
			int maxIterationsToRun = 30;
			int percentageOfStableNodes = 99;
			String partitionFile = "";
			boolean partitionIsFromFile = false;
			
			if (args.length > 0)
				pathToGraph = args[0];		
			
			if (args.length > 1)
				outputPath = args[1];		
			
			if (args.length > 2)
				betas = ParseDoubleArray(args[2]);		
				
			if (args.length > 3)
				 alpha = Double.parseDouble(args[3]);
			
			if ( args.length > 4)
				 iteratioNumToStartMerge = Integer.parseInt(args[4]);
			
			if ( args.length > 5)
				maxIterationsToRun = Integer.parseInt(args[5]);
			
			if(args.length > 6){
				percentageOfStableNodes = Integer.parseInt(args[6]);
				if(percentageOfStableNodes<1 || percentageOfStableNodes>100){
					throw(new RuntimeException("param at location 6 is percentageOfStableNodes. You gave: " + percentageOfStableNodes +"  which is not <1 or >100."));
				}
			}
			
			if(args.length > 7){
				partitionFile = args[7];
				partitionIsFromFile = true;
			}
			
			String betasString = "";
			for (double d: betas){
				betasString = betasString + d + " ";
			}
			
			System.out.println("pathToGraph:             "+pathToGraph);      
			System.out.println("outputPath:              "+outputPath);
			System.out.println("betas:                   "+betasString);
			System.out.println("alpha:                   "+alpha);
			System.out.println("iteratioNumToStartMerge: "+iteratioNumToStartMerge);
			System.out.println("maxIterationsToRun:      "+maxIterationsToRun);
			System.out.println("percentageOfStableNodes: "+percentageOfStableNodes);
			if(partitionIsFromFile){
				System.out.println("partitionFile:           "+ partitionFile);
			}
			System.out.println("");
			
			OSCDWithlevels oscd = null;
			if (partitionIsFromFile){
				oscd = new OSCDWithlevels(pathToGraph, partitionFile, betas,alpha,outputPath, iteratioNumToStartMerge, maxIterationsToRun, percentageOfStableNodes);
			}
			else{
				oscd = new OSCDWithlevels(pathToGraph,betas,alpha,outputPath, iteratioNumToStartMerge, maxIterationsToRun,percentageOfStableNodes);
			}
			oscd.FindCommunities();
		}
	}

	private static double[] ParseDoubleArray(String string) {
		String[] parts = string.split(",");
		double[] ans= new double[parts.length];
	    int i=0;
	    for(String str:parts){
	    	ans[i]=Double.parseDouble(str);
	        i++;
	    }
		return ans;
	}
}