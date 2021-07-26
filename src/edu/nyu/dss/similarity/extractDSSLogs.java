package edu.nyu.dss.similarity;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import au.edu.rmit.trajectory.clustering.kpaths.Util;

public class extractDSSLogs {
	static String folderString = "./logs/spadas/efficiency/";
	
	
	public static void main(String[] args) throws FileNotFoundException, IOException {
		extractLogs(5);
	}
	
	
	static // extract the data from all log files for the end-to-end processing
	void extractLogs(int N) throws FileNotFoundException, IOException {
		File folder = new File(folderString);
		for (File f : folder.listFiles()) {
		    if (f.getName().endsWith(".csv")) {
		        f.delete();
		    }
		}
		for(int i=1; i<=N; i++) {
			String logs = readAverage(folderString+i+"-r-topk-grid.txt");
			Util.rewrite(folderString+"r-topk-grid.csv", logs+"\n");
			
			
			logs = readAverage(folderString+i+"-r-haus-appro.txt");
			Util.rewrite(folderString+"r-topk-haus.csv", logs+"\n");
			
			
			logs = readAverage(folderString+i+"-r-range.txt");
			Util.rewrite(folderString+"r-range.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-r-range-deep.txt");
			Util.rewrite(folderString+"r-range-deep.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-N-range.txt");
			Util.rewrite(folderString+"N-range.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-N-topk-haus.txt");
			Util.rewrite(folderString+"N-topk-haus.csv", logs+"\n");
			
			
			
		//	logs = readAverage(folderString+i+"-d-topk-haus.txt");
		//	Util.write(folderString+"d-topk-haus.txt", logs+"\n");
			
			
			
		//	logs = readAverage(folderString+i+"-d-range.txt");
		//	Util.write(folderString+"d-range.txt", logs+"\n");
			
		//	logs = readAverage(folderString+i+"-d-range-deep.txt");
		//	Util.write(folderString+"d-range-deep.txt", logs+"\n");
			
			
			logs = readAverage(folderString+i+"-k-range.txt");
			Util.rewrite(folderString+"k-range.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-k-topk-haus.txt");
			Util.rewrite(folderString+"k-topk-haus.csv", logs+"\n");
			
		//	logs = readAverage(folderString+i+"-k-range-deep.txt");
		//	Util.rewrite(folderString+"k-range-deep.txt", logs+"\n");
			
			logs = readAverage(folderString+i+"-n-haus.txt");
			Util.rewrite(folderString+"n-haus.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-n-join-outlier.txt");
			Util.rewrite(folderString+"n-join-outlier.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-f-haus.txt");
			Util.rewrite(folderString+"f-haus.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-f-range.txt");
			Util.rewrite(folderString+"f-range.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-f-topk.txt");
			Util.rewrite(folderString+"f-topk.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-N-MemorySize.txt");
			Util.rewrite(folderString+"N-MemorySize.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-N-ConstructionCost.txt");
			Util.rewrite(folderString+"N-ConstructionCost.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-M-N-topk.txt");
			Util.rewrite(folderString+"M-N-topk.csv", logs+"\n");
			
			logs = readAverage(folderString+i+"-M-N-range.txt");
			Util.rewrite(folderString+"M-N-range.csv", logs+"\n");
		}
		
		String logs = readAverage(folderString+"5-d-topk-haus.txt");
		Util.rewrite(folderString+"d-topk-haus.csv", logs+"\n");
		
		logs = readAverage(folderString+"5-d-range.txt");
		Util.rewrite(folderString+"d-range.csv", logs+"\n");
		
		logs = readAverage(folderString+"5-d-haus.txt");
		Util.rewrite(folderString+"d-haus.csv", logs+"\n");
		
		logs = readAverage(folderString+"5-M-d-topk.txt");
		Util.rewrite(folderString+"M-d-topk.csv", logs+"\n");
		
		logs = readAverage(folderString+"5-M-d-range.txt");
		Util.rewrite(folderString+"M-d-range.csv", logs+"\n");
		
		logs = readAverage(folderString+"5-d-MemorySize.txt");
		Util.rewrite(folderString+"d-MemorySize.csv", logs+"\n");
		
		logs = readAverage(folderString+"5-d-ConstructionCost.txt");
		Util.rewrite(folderString+"d-ConstructionCost.csv", logs+"\n");
		
		extractIndexLogs();
	}
	
	/*
	 * combining each dataset's log for plotting
	 */
	static String readAverage(String filename) throws FileNotFoundException, IOException{
		File file = new File(filename);
		double[] logs = null;
		int count = 0;
		int length = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			while ((strLine = br.readLine()) != null) {
				String[] splitString = strLine.split(",");
				length = splitString.length;
				if(logs==null)
					logs = new double[splitString.length];
				for(int i=0; i<splitString.length; i++)
					logs[i] += Double.valueOf(splitString[i]);
				count++;
			}
		}
		String content = "";
		for(int i=0; i<length; i++)
			content += logs[i]/count+",";
		return content;
	}
	
	static double[][] readArrays(String filename) throws FileNotFoundException, IOException{
		File file = new File(filename);
		double measures[][] = new double[5][4*5+1];
		int count = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String strLine;
			while ((strLine = br.readLine()) != null) {
				String[] splitString = strLine.split(",");
			//	System.out.println(strLine);
				if(count>=5)
					break;
				for(int i=0; i<splitString.length; i++)
					measures[count][i] = Double.valueOf(splitString[i]);
				count++;
			}
		}
		return measures;
	}
	
	/*
	 * read the index logs to extract infor about index construction
	 */
	static void extractIndexLogs() throws FileNotFoundException, IOException {
		// extract the logs from generated files, and choose the first value 
		double[][] index = readArrays(folderString+"N-ConstructionCost.csv");
		double[][] topkrange = readArrays(folderString+"r-range.csv");
		double[][] topkgbo = readArrays(folderString+"r-topk-grid.csv");
		double[][] topkhaus = readArrays(folderString+"r-topk-haus.csv");
		double[][] join = readArrays(folderString+"n-join-outlier.csv");
		double[][] union = readArrays(folderString+"r-range-deep.csv");
		
		for(int i=0; i<5; i++) {
			String content = "";
			content += Integer.toString(i+1)+",";
			content += index[i][4]+",";
			content += topkrange[i][4]+",";
			content += topkgbo[i][4]+",";
			content += topkhaus[i][4]+",";
			content += join[i][4]+",";
			content += union[i][4]+",";
			Util.rewrite(folderString+"overall.csv", content+"\n");
		}
	}
}
