package edu.nyu.dss.similarity;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.PriorityQueue;

import org.apache.commons.io.FileUtils;

import au.edu.rmit.trajectory.clustering.kmeans.indexAlgorithm;
import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import au.edu.rmit.trajectory.clustering.kpaths.Util;


public class EffectivenessStudy {
	static boolean dimSelected[];// indicate which dimension is not selected, for inferring the radius
	static boolean dimensionAll = true;// indicate that all dimensions are selected.
	
	/*
	 * combine two integers to produce a new value
	 */
	public static int combine(int aid, int bid, int lengtho){
		int length = lengtho;
		int[] a =new int[length];
		int[] b =new int[length];
		while(length-- >= 1){
			a[length] = aid%2;
			aid /=2;
			b[length] = bid%2;
			bid /=2;
		}
		int com[] = new int[2*lengtho];
		for(int i = 0; i<lengtho; i++){
			com[2*i]= a[i];
			com[2*i+1] = b[i];
		}
		return bitToint(com, 2*lengtho);
	}
	
	/*
	 * generate the z-curve code
	 */
	public static int bitToint(int[] a, int length){
		int sum = 0;
		for(int i=0; i<length; i++){
			sum += a[i]*Math.pow(2, length-i-1);
		}
		return sum;
	}
	
	
	public static void SerializedZcurve(String file, HashMap<Integer, ArrayList<Integer>> result) {
		try {
			FileOutputStream fos = new FileOutputStream(file);
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(result);
			oos.close();
			fos.close();
			System.out.println("Serialized result HashMap data is saved in hashmap.ser");
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	/*
	 * load the zcurve into memory
	 */
	static HashMap<Integer, ArrayList<Integer>> deSerializationZcurve(String file) {
		HashMap<Integer, ArrayList<Integer>> result;
		try {
			FileInputStream fis = new FileInputStream(file);
			ObjectInputStream ois = new ObjectInputStream(fis);
			result = (HashMap) ois.readObject();
			ois.close();
			fis.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
			return null;
		} catch (ClassNotFoundException c) {
			System.out.println("Class not found");
			c.printStackTrace();
			return null;
		}
		return result;
	}
	
	/*
	 * resolution is the parameter determine how many cells will be generated, only for two dimension here
	 */
	public static HashMap<Integer, ArrayList<Integer>> storeZcurve(double minx, double miny, double range, 
			int dimension, Map<Integer, double[][]> datamap, int resolution, String zcodeFile) {
		HashMap<Integer, ArrayList<Integer>> zcodeMap = new HashMap<Integer, ArrayList<Integer>>();// using the z-curve to store the encoding based on binary of x, y
		int numberCells = (int) Math.pow(2, resolution);
		double unit = range/numberCells;
		for(int datasetid: datamap.keySet()) {
			ArrayList<Integer> zcodeaArrayList = new ArrayList<Integer>();
			double [][]dataset = datamap.get(datasetid);
			for(int i=0; i<dataset.length; i++) {
				int x = (int)((dataset[i][0]-minx)/unit);
				int y = (int)((dataset[i][1]-miny)/unit);
				int zcode = combine(x,y,resolution);
			//	System.out.println(zcode);
				if(!zcodeaArrayList.contains(zcode))
					zcodeaArrayList.add(zcode);
			}
			zcodeMap.put(datasetid, zcodeaArrayList);
		}
		SerializedZcurve(zcodeFile, zcodeMap);
		return zcodeMap;
	}
	
	// to do, 
	void ComputeEMD(ArrayList<Integer> a, ArrayList<Integer> b, double unit, int resolution) {
		// 1. exact,
		
		// 2. approximate, Table
		
	}
	
	void LowerBoundEMD() {
		// compute the lower bound
		
		// Table 1.
		
	}
	
	void BuildIndex(indexNode root) {
		// compute pivot and radius
		
		// 
	}
	
	void TopKEMD() {
		// 1. linear scanning, compute the lower bound first, then the exact, or approximate
		
		// 2. use index to prune, 
	}
	
	public static HashMap<Integer, ArrayList<Integer>> storeZcurveFile(double minx, double miny, double range, 
			int dimension, int resolution, String zcodeFile, int limit, Map<Integer, String> DatasetMapping) throws FileNotFoundException, IOException {
		HashMap<Integer, ArrayList<Integer>> zcodeMap = new HashMap<Integer, ArrayList<Integer>>();// using the z-curve to store the encoding based on binary of x, y
		int numberCells = (int) Math.pow(2, resolution);
		double unit = range/numberCells;
		for(int datasetid = 1; datasetid<=limit; datasetid++) {
			ArrayList<Integer> zcodeaArrayList = new ArrayList<Integer>();
			double [][]dataset;
			dataset = Framework.readSingleFile(DatasetMapping.get(datasetid));
			for(int i=0; i<dataset.length; i++) {
				int x = (int)((dataset[i][0]-minx)/unit);
				int y = (int)((dataset[i][1]-miny)/unit);
				int zcode = combine(x,y,resolution);
				if(!zcodeaArrayList.contains(zcode))
					zcodeaArrayList.add(zcode);
			}
			zcodeMap.put(datasetid, zcodeaArrayList);
		}
		SerializedZcurve(zcodeFile, zcodeMap);
		return zcodeMap;
	}
	
	/*
	 * scan every dataset to find the top-k
	 */
	static HashMap<Integer, Double> topkAlgorithmZcurveHashMap(Map<Integer, ArrayList<Integer>> bitMap, ArrayList<Integer> query, int k) {
		HashMap<Integer, Double> result = new HashMap<Integer, Double>();
		HashMap<Integer, PriorityQueue<queueMain>> queues = new HashMap<Integer, PriorityQueue<queueMain>>();
		for(int datasetid: bitMap.keySet()) {
			ArrayList<Integer> arrayList = bitMap.get(datasetid);
			int counter = 0;
			for(int queryz: query) {
				if(arrayList.contains(queryz)) {
					counter++;
				}
			}
			double reverseoverlapRatio = 1 - counter/(double)query.size();
			result = Search.holdingTopK(result, datasetid, reverseoverlapRatio, k, queues, null);
		}
		return result;
	}
	
	// store the top-k result into a file for further operation
	public static void SerializedResultTable(String file, HashMap<Integer, Double> result) {
		try {
			FileOutputStream fos = new FileOutputStream(file);
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(result);
			oos.close();
			fos.close();
			System.out.println("Serialized result HashMap data is saved in hashmap.ser");
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	static HashMap<Integer, Double> deSerializationResult(String file) {
		HashMap<Integer, Double> result;
		try {
			FileInputStream fis = new FileInputStream(file);
			ObjectInputStream ois = new ObjectInputStream(fis);
			result = (HashMap) ois.readObject();
			ois.close();
			fis.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
			return null;
		} catch (ClassNotFoundException c) {
			System.out.println("Class not found");
			c.printStackTrace();
			return null;
		}
		return result;
	}
	
	/*
	 * mapping table 
	 */
	static Map<Integer, String> deSerialization(String file) {
		Map<Integer, String> argoDataMap;
		try {
			FileInputStream fis = new FileInputStream(file);
			ObjectInputStream ois = new ObjectInputStream(fis);
			argoDataMap = (HashMap) ois.readObject();
			ois.close();
			fis.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
			return null;
		} catch (ClassNotFoundException c) {
			System.out.println("Class not found");
			c.printStackTrace();
			return null;
		}
		return argoDataMap;
	}
	
	// join: we increase the threshold to observe the matching ratio, in spatial join-based augmentation
	static void writeJoinRatio(double joinThreshold[], ArrayList<Double> distancearray, String file) {
		File myObj = new File(file); 
    	myObj.delete();
		for(double threshold: joinThreshold) {
			int stoppoint = 0;
			for(double distance: distancearray) {
				stoppoint++;
				if(distance<threshold) {
					double joinratio = (distancearray.size()-stoppoint)/(double)distancearray.size();
					Util.write(file, joinratio+"\n");
					break;
				}
			}
		}
	}
	
	/*
	 * split the Chicago dataset by 15mins
	 */
	static void splitChicagoTrips(long limit) throws FileNotFoundException, IOException {
		String filename = "/Volumes/SSD/Taxi_Trips_Chicago.csv";
		File folder = new File(filename);
		try (BufferedReader br = new BufferedReader(new FileReader(folder))) {
			String strLine;
			int counter = 0;
			while ((strLine = br.readLine()) != null) {			
				counter++;
			//	if(counter<53328075)
			//		continue;
			//	System.out.println(strLine);
				String[] splitString = strLine.split(",");
				if(counter==1)
					continue;
				if(counter%10000==0)
					System.out.println(counter);
				// extract time, and put in the same file that share the same hour,
				String starttime =  splitString[2];
				if(splitString.length<=21)
					continue;
				String starttimes[] =  starttime.split(" ");
				String starttimess[] =  starttimes[0].split("/");
				String starttimess1[] =  starttimes[1].split(":");
				int hour = Integer.valueOf(starttimess1[0]);
				if(starttimes[2].equals("PM")) {
					hour += 12;
				}
				
				String tripfile = starttimess[2]+starttimess[0]+starttimess[1]+hour+starttimess1[1];
				String tripsecond = splitString[4];
				String tripmile = splitString[5];
				String fare = splitString[10];
				String tips = splitString[11];
				String tolls = splitString[12];
				String Extras = splitString[13];
				String triptotal = splitString[14];
				String pickuplatitude = splitString[17];
				String pickupllongtigutde = splitString[18];
					
				String droplatitude = splitString[20];
				String dropllongtigutde = splitString[21];
				String cleandataString = tripsecond+","+tripmile+","+fare+","+tips+","+tolls+","+Extras+","+triptotal+","+pickuplatitude+","+pickupllongtigutde+
						","+droplatitude+","+dropllongtigutde;
			//	System.out.println(cleandataString);
				
				Util.write("/Users/sw160/Desktop/chicago15mins/"+tripfile+".txt", cleandataString+"\n");//write into the same file that share the same pickup time,
			}
		}
	}

	
	/*
	 * this is for interactive search, where index of data lake and each dataset are restored, instead of building from scratch.
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {
	//	splitChicagoTrips(10000000);
	//	combineAllFiles("/Volumes/SSD/chicago15mins39020000/", "/Volumes/SSD/chicago15mins/"); //combine the new files
			
		/*
		 * test the index of whole data lake
		 */
		double joinThreshold[] = {50, 100, 150, 200, 250};
		indexAlgorithm<Object> indexDSS = new indexAlgorithm<>();
		System.out.println("indexing datasets");
		Map<Integer, indexNode> datalakeIndex = null;
		indexDSS.setGloabalid();
		String indexString = "./index/dss/index/argo/";
		int dimension = 3;
		dimSelected = new boolean[dimension];
		dimSelected[2] = true;
		dimensionAll = false;
		Map<Integer, String> argoDataMap = deSerialization(indexString+"mapping.ser");
		System.out.print(argoDataMap.size());
<<<<<<< HEAD
//		modify the foldername to datalake1-10-3
//		datalakeIndex = indexDSS.restoreDatalakeIndex(indexString+"datalake.txt", dimension);
		datalakeIndex = indexDSS.restoreDatalakeIndex(indexString+"datalake1-10-3.txt", dimension);
=======
		datalakeIndex = indexDSS.restoreDatalakeIndex(indexString+"datalake.txt", dimension);
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
		//create signature, bounding box after the restoring
		Map<Integer, Map<Integer, indexNode>> datasetIndex = null;
//		datasetIndex = indexDSS.restoreIndex(indexString, dimension, null); //restore all the index, we do not need this when we have data lake index
		Map<Integer, double[][]> dataMapPorto = null;
		
		int queryID = Integer.valueOf(args[0]);
		int k = Integer.valueOf(args[1]);
		int joinid = Integer.valueOf(args[2]);
		int topkAlgorithmOption = Integer.valueOf(args[5]);
		// we can also specify one file in a query folder
		HashMap<Integer, Double> result = new HashMap<Integer, Double>();
		Map<Integer, indexNode> queryindexMap, dataindexMap;
		System.out.println(argoDataMap.get(queryID));
		if(datasetIndex==null) {
			queryindexMap = indexDSS.restoreSingleIndex(indexString, queryID, dimension);
			// or recreate the index for each single
		}else {
			queryindexMap = datasetIndex.get(queryID);
		}
		Map<Integer, double[][]> dataMap = new HashMap<Integer, double[][]>();
		if(!args[4].equals("")) {
			// or read a dataset from folder, and search, then put the related dataset to 
		}
		String zcodeSer= "./logs/spadas/"+args[4]+"/argo_bitmaps.ser";
		File tempFile = new File(zcodeSer);
		HashMap<Integer, ArrayList<Integer>> zcodemap = new HashMap<Integer, ArrayList<Integer>>();
		AdvancedHausdorff.setWeight(3, indexString+"weight.ser");//set weight on all dimensions, for point distance computation
		if(!tempFile.exists()) {
			// read the porto
			double minx = 0, miny=500, range = 4600;
			int resolution = 8;
			String aString = "/Users/sw160/Desktop/argoverse-api/dataset/train/data";
<<<<<<< HEAD
//			String aString = "./dataset/argoverse/opendd_data";
=======
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
			int limit = 10000;
			File folder = new File(aString);
			dataMap = Framework.readFolder(folder, limit);
			zcodemap = storeZcurve(minx, miny, range, 2, dataMap, resolution, zcodeSer);
		}else {
			zcodemap = deSerializationZcurve(zcodeSer);
		}
		//create signature and max covered points on the restored index
		datalakeIndex.get(1).setMaxCovertpoint(datalakeIndex);//
		datalakeIndex.get(1).buildsignarture(zcodemap,datalakeIndex);//
		System.out.println("max:"+datalakeIndex.get(1).getmaxCoverpoints());
	//	System.out.println("query is "+queryindexMap.get(1).getMBRmax()[0]+","+queryindexMap.get(1).getMBRmax()[1]+","+queryindexMap.get(1).getMBRmax()[2]);
	//	System.out.println("query is "+queryindexMap.get(1).getMBRmin()[0]+","+queryindexMap.get(1).getMBRmin()[1]+","+queryindexMap.get(1).getMBRmin()[2]);
		int capacity = 30;
		double weight[] = null;
		if (args[3].equals("topk")) {
			// read the query datasets, read the dataset, if both are not avaiable, we read.
			long startTime1 = System.nanoTime();
			if(topkAlgorithmOption==1)//hausdorff
				result = Search.pruneByIndex(dataMapPorto, null, queryindexMap.get(1), queryID, dimension, null,
					datasetIndex, queryindexMap, datalakeIndex, argoDataMap, k, indexString, dimSelected, dimensionAll, 0, capacity, weight, false, null);
			else if(topkAlgorithmOption==2)//grid-based overlap, scan every datasets
				result = topkAlgorithmZcurveHashMap(zcodemap, zcodemap.get(queryID), k);
			else if(topkAlgorithmOption==3)//rank by number of covered points
				result = Search.rangeQueryRankingNumberPoints(datalakeIndex.get(1), result, queryindexMap.get(1).getMBRmax(), 
						queryindexMap.get(1).getMBRmin(), Double.MAX_VALUE, k, null, dimension, dataMapPorto, datalakeIndex, 
						datasetIndex, argoDataMap, indexString, indexDSS, dimSelected, dimensionAll, null, null, capacity, weight, false);
			else if(topkAlgorithmOption==4)//rank by area
				result = Search.rangeQueryRankingArea(datalakeIndex.get(1), result, queryindexMap.get(1).getMBRmax(), 
						queryindexMap.get(1).getMBRmin(), Double.MAX_VALUE, k, null, dimension, datalakeIndex, dimSelected, dimensionAll);
			else if(topkAlgorithmOption==5) {//rank by grid-based overlap, fast
				int []queryzcurve = new int[zcodemap.get(queryID).size()];
				for(int i=0; i<zcodemap.get(queryID).size(); i++)
					queryzcurve[i] = zcodemap.get(queryID).get(i);
				result = Search.gridOverlap(datalakeIndex.get(1), result, queryzcurve, Double.MAX_VALUE, k, null, datalakeIndex);
			}
			// call the overlap ratio
			long endtime = System.nanoTime();
			System.out.println("top-1 search costs: " + (endtime - startTime1) / 1000000000.0);
			System.out.println(result);
			SerializedResultTable("./logs/spadas/"+args[4]+"/result.ser", result);
		} else if(args[3].equals("join")){//conduct join
			result = deSerializationResult("./logs/spadas/"+args[4]+"/result.ser");
			int resultid = 0;//result.entrySet().iterator().next().getKey();// we select the k-th results.
			int counter = 10;
			for (int id : result.keySet()) {
				counter--;
				if (counter < joinid) {
					resultid = id;
					break;
				}
			}
			System.out.println(result);
			System.out.println(resultid);
			if (datasetIndex == null) {
				dataindexMap = indexDSS.restoreSingleIndex(indexString, resultid, dimension);
			} else {
				dataindexMap = datasetIndex.get(resultid);
			}
			// conduct join afterward, for effectiveness
			PriorityQueue<queueMain> aHeaps = new PriorityQueue<queueMain>();
			PriorityQueue<queueSecond> secondHeaps = new PriorityQueue<queueSecond>();
			queueSecond secondq = new queueSecond(dataindexMap.get(1), 0, 0);
			secondHeaps.add(secondq);
			queueMain mainq = new queueMain(queryindexMap.get(1), secondHeaps, Double.MAX_VALUE, 0);
			aHeaps.add(mainq);
			double[][] querydata, dataset;
			querydata = Framework.readSingleFile(argoDataMap.get(queryID));// read from files
			dataset = Framework.readSingleFile(argoDataMap.get(resultid));
			ArrayList<Double> distancearray = Join.IncrementalJoin(querydata, dataset, dimension, queryindexMap.get(1), dataindexMap.get(1), 0, 0, 0.01,
					false, 0, false, Double.MAX_VALUE, aHeaps, queryindexMap, dataindexMap, args[4], dimSelected, dimensionAll);
			writeJoinRatio(joinThreshold, distancearray, "./logs/spadas/"+args[4]+"/joinRatio.txt");//evaluate the matching ratio, 
		}else if(args[3].equals("union")) {
			// union operation, put the dataset into a single folder for training, and prediction argoverse
			result = deSerializationResult("./logs/spadas/"+args[4]+"/result.ser");
			FileUtils.cleanDirectory(new File(args[6])); //clear the folder
			for (int id : result.keySet()) {				
				Framework.UnionFile(argoDataMap.get(id), args[6]);
			}
			File myObj = new File("./logs/spadas/"+args[4]+"/.DS_Store"); 
	    	myObj.delete();
		}else if(args[3].equals("union-range")) {
			// if the query is a range, we return all the points in the region of the selected datasets
			// 
		}
	}	
}
