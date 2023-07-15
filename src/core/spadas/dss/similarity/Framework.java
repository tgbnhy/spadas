package core.spadas.dss.similarity;

import java.io.*;
import java.nio.file.Files;

import java.util.*;
import java.util.stream.Stream;

import core.spadas.dss.similarity.onethreeseven.Kneedle;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeMap;


import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.tuple.Pair;

import tree.trajectory.clustering.kmeans.indexAlgorithm;
import tree.trajectory.clustering.kmeans.indexNode;
import tree.trajectory.clustering.kpaths.Util;


public class Framework {
    /*
     * data set
     */
    //	static int dimension = 3;
//	originally 10
    /*
     * argo-data: 30
     * argo-poi: 10
     * beijing: 100
     * porto: 150
     * shapenet: 50
     * nyu: 2000
     * */
    static int capacity = 100;
    static String edgeString, nodeString, distanceString, fetureString, indexString = "";
    static String zcodeSer;
    static String aString = "/Users/sw160/Desktop/argoverse-api/dataset/train/data";
    static Map<Integer, String> datasetIdMapping = new HashMap<Integer, String>();//integer
    static Map<Integer, double[][]> dataMapPorto = new HashMap<Integer, double[][]>();
    static TreeMap<Integer, Integer> countHistogram = new TreeMap<Integer, Integer>();

    static TreeMap<Integer, Integer> countHistogramNew = new TreeMap<Integer, Integer>();
    static int datalakeID;// the lake id
    static String folderString = ".";
    //    一个存放数据集集中每个数据集包含数据点数目的列表，用来绘制直方图，和前面的countHistogram作用类似，但适用于nyu和btaa数据集
//    虽然需要排序但每次由于调用次数不多所以就随便排不优化了
    static List<Integer> pointCntList = new ArrayList<>();

    //	record covered points number, <node's number, covered points number>
    static Map<Integer, Integer> coveredPointsMap = new TreeMap<>();
    //	record radius, <default index of leaf indexNode, radius>
    static List<Double> radiusList = new ArrayList<>();
    static PriorityQueue<Double> radiusQueue = new PriorityQueue<>();

    /*
     * z-curve for grid-based overlap
     */
    static double minx = 0, miny = 500;// the range of dataset space, default values for Argo dataset
    static double spaceRange = 4600;
    static int resolution = 3; // also a parameter to test the grid-based overlap, and the approximate hausdorff, and range query
    static double error = 0.0;

    /*
     * index
     */
    static indexAlgorithm<Object> indexDSS = new indexAlgorithm<>();
    static Map<Integer, indexNode> indexMap;// root node of dataset's index
    static Map<Integer, indexNode> datalakeIndex = null;// the global datalake index
    static ArrayList<indexNode> indexNodes = new ArrayList<indexNode>(); // store the root nodes of all datasets in the lake
    static ArrayList<indexNode> indexNodesAll;
    static indexNode datasetRoot = null; // the root node of datalake index in memory mode
    static Map<Integer, Map<Integer, indexNode>> datasetIndex; // restored dataset index
    static HashMap<Integer, ArrayList<Integer>> zcodemap;
    static Map<Integer, Pair<Double, PriorityQueue<queueMain>>> resultPair =
            new HashMap<Integer, Pair<Double, PriorityQueue<queueMain>>>();
    /*
     * dataset features
     */
    static Map<Integer, double[]> featureMap = new HashMap<Integer, double[]>(); //storing the features for each map

    /*
     * parameters
     */
    static double weight[] = null; // the weight in all dimensions
    static int dimension = 2;

    static int limit = 100000; //number of dataset to be searched
    static int NumberQueryDataset = 1; // number of query datasets
    static boolean dimNonSelected[];// indicate which dimension is not selected
    static boolean dimensionAll = true;// indicate that all dimensions are selected.

    static int ks[] = {10, 20, 30, 40, 50};
    static int ds[] = {2, 4, 6, 8, 10};
    static double ns[] = {0.0001, 0.001, 0.01, 0.05, 0.1};// number of datasets indexed
    static int ss[] = {1, 10, 100, 500, 1000}; //number of datasets to be combined together, i.e., the scale of dataset
    static int rs[] = {3, 4, 5, 6, 7};// resolution to verify the range query, grid overlap, approximate hausdorff
    static int fs[] = {10, 20, 30, 40, 50};


    //	max radius
    static double maxRadius = 0;
    //	number of leaf node whose radius is within max radius
    static int normalNodes = 0;
    //	test time
    static double buildIndexTreeTime = 0;
    static double deleteOutliersTime = 0;
    static int pointsRecorded = -1;

    /*
     * settings
     */
    static boolean storeIndexMemory = true;// check whether need to use the index in memeory or load from disk

    //	change false to true
//	change true to false
    static boolean saveDatasetIndex = false;// whether to store each dataset index as file
    //	change false to true
    static boolean storeAllDatasetMemory = true; // store all the datasets in the datasetPorto

    static boolean zcurveExist = false;
    static boolean buildOnlyRoots = false; // for building fast datalake index
    static boolean testFullDataset = true; //

    /*
     * information to compare all the algorithms
     */
    static int numberofComparisons = 4;// all the comparisons
    static int numberofParameterChoices = 5; //
    static double timeMultiple[][];// the overall running time
    static double indexMultiple[][];
    static int fileNo = 1;

    /*
     * read the shapenet dataset
     */
    static Map<Integer, double[][]> readShapeNet(File folder, Map<Integer, double[][]> datasetMap, int limit) throws IOException {
        File[] fileNames = folder.listFiles();
        fileNo = 1;
        datasetMap = null;
        if (storeAllDatasetMemory)
            datasetMap = new HashMap<Integer, double[][]>();
        datasetIdMapping = new HashMap<Integer, String>();
        for (File file : fileNames) {
            if (file.isDirectory()) {
                //read the files inside
                readShapeNet(file, datasetMap, limit);
            } else {
                if (!file.getName().contains("_") && !file.getName().contains("name")) {
                    System.out.println(file.getName());
                    //	System.out.println(file.getName());
                    datasetIdMapping.put(fileNo, file.getName());
                    long lineNumber = 0;
                    try (Stream<String> lines = Files.lines(file.toPath())) {
                        lineNumber = lines.count();
                    }
                    double[][] a = new double[(int) lineNumber][];
                    int i = 0;
                    try (BufferedReader br = new BufferedReader(new FileReader(file))) {
                        String strLine;
                        String[] splitString;
                        while ((strLine = br.readLine()) != null) {
                            // System.out.println(strLine);
                            splitString = strLine.split(" ");
                            a[i] = new double[3];
                            a[i][0] = new Double(splitString[0]);
                            a[i][1] = new Double(splitString[0]);
                            a[i][2] = new Double(splitString[0]);
//                            a[i][0] = Double.parseDouble(splitString[0]);
//                            a[i][1] = Double.parseDouble(splitString[1]);
//                            a[i][2] = Double.parseDouble(splitString[2]);
                            i++;
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    //	System.out.println(i+","+lineNumber);

                    if (countHistogram.containsKey(a.length))// count the histograms
                        countHistogram.put(a.length, countHistogram.get(a.length) + 1);
                    else
                        countHistogram.put(a.length, 1);
                    if (storeAllDatasetMemory == true) {
                        datasetMap.put(fileNo, a);
                    }
                    long startTime = System.nanoTime();
                    if (storeIndexMemory) {
                        createDatasetIndex(fileNo, a);
                    }
                    long endTime = System.nanoTime();
                    buildIndexTreeTime += (endTime - startTime) / 1000000000.0;
                    if (!zcurveExist)
                        storeZcurve(a, fileNo);
                    //	Path sourceDirectory = file.toPath();// Paths.get("/Users/personal/tutorials/source");
                    //  Path targetDirectory = Paths.get("/Users/sw160/Desktop/spadas-dataset/Shapenet-allInOne/"+fileNo+".txt");
                    //  write("/Users/sw160/Desktop/spadas-dataset/Shapenet-allInOne/namemapping.txt", fileNo+","+file.getName()+"\n");
                    //  Files.copy(sourceDirectory, targetDirectory); // write into the same folder
                    fileNo++;
                    if (fileNo > limit)
                        break;

                }
            }
        }
        if (fileNo < limit)
            limit = fileNo;
        return datasetMap;
    }

    /*
     * read the trajectory database, e.g., Porto, T-drive
     */
    static Map<Integer, double[][]> readTrajectoryData(File file, int limit) throws IOException {
        Map<Integer, double[][]> datasetMap = null;
        if (storeAllDatasetMemory)
            datasetMap = new HashMap<Integer, double[][]>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String strLine;
            int line = 1;
            while ((strLine = br.readLine()) != null) {
                String[] splitString = strLine.split(",");
                double[][] a = new double[splitString.length / 2][];
                for (int i = 0; i < splitString.length; i++) {
                    if ((i + 1) % 2 == 0) {
                        a[i / 2] = new double[2];
                        for (int j = i - 1; j < i + 1; j++) {
                            a[i / 2][j - i + 1] = Double.valueOf(splitString[j]);
                        }
                    }
                }
                if (countHistogram.containsKey(splitString.length / 2))
                    countHistogram.put(splitString.length / 2, countHistogram.get(splitString.length / 2) + 1);
                else
                    countHistogram.put(splitString.length / 2, 1);
                if (storeAllDatasetMemory)
                    datasetMap.put(line, a);
                long startTime = System.nanoTime();
                if (storeIndexMemory) {
                    createDatasetIndex(line, a);
                }
                long endTime = System.nanoTime();
                buildIndexTreeTime += (endTime - startTime) / 1000000000.0;
                if (!zcurveExist)
                    storeZcurve(a, line);
                line++;
                if (line > limit + NumberQueryDataset)
                    break;

            }
            if (line < limit)
                limit = line;

        }
        return datasetMap;
    }

    static Map<Integer, double[][]> readNyu(File folder, int limit) {
        File[] fileNames = folder.listFiles();
//        int fileNo = 1;
//        if (storeAllDatasetMemory)
//            dataMapPorto = new HashMap<Integer, double[][]>();
        for (File file : fileNames) {
            if (file.isDirectory()) {
                readNyu(file, limit);
            } else {
                try {
                    String a = file.getName();
                    datasetIdMapping.put(fileNo, a);
//                    根据数据集集的名称是否包含中文来分类读取
                    char firstChar = a.charAt(0);
                    if (Character.isIdeographic((int) firstChar)) {
                        readContentBaidu(file);
                    } else {
                        readContentNyu(file);
                    }
                    fileNo++;
//                    Collections.sort(pointCntList);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            if (fileNo > limit)
                break;
        }
        return dataMapPorto;
    }

    public static void readContentNyu(File file) throws IOException {
        System.out.println("read file " + file.getCanonicalPath());
        long lineNumber = 0;
        try (Stream<String> lines = Files.lines(file.toPath())) {
            lineNumber = lines.count();
        }
        List<double[]> l = new ArrayList<>();
//        double[][] a = new double[(int) lineNumber - 1][];// a little different
        int i = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String strLine;
            while ((strLine = br.readLine()) != null) {
                // System.out.println(strLine);
                String[] splitString = strLine.split(",");
//				3 for folder "data"
//				0 for folder "poi"
                String aString = splitString[0];
                if (aString.matches("-?\\d+(\\.\\d+)?")) {// only has float
//                    double[] b = new double[2];
//                    a[i] = new double[2];
//                    a[i][0] = Double.parseDouble(splitString[0]);
//                    a[i][1] = Double.parseDouble(splitString[1]);
                    double[] b = new double[2];
                    b[0] = new Double(splitString[0]);
//                    b[0] = Double.parseDouble(splitString[0]);
//                    b[1] = Double.parseDouble(splitString[1]);
                    b[1] = new Double(splitString[1]);
//                    if (b[0] < -180 || b[0] > 180 || Math.abs(b[0]) < 0.001 || b[1] < -90 || b[1] > 90 || Math.abs(b[1]) < 0.001) {
//                        continue;
//                    }
                    if (Math.abs(b[0]) < 0.01 || Math.abs(b[1]) < 0.01) {
                        continue;
                    }
//					a[i][2] = Double.valueOf(splitString[0]);
//					a[i][0] = Double.valueOf(splitString[0]);
//					a[i][1] = Double.valueOf(splitString[1]);
                    l.add(b);
                    i++;
                }
            }
//			long endTime = System.nanoTime();
//			buildIndexTreeTime += (endTime - startTime) / 1000000000.0;
        }
        double[][] a = l.toArray(new double[][]{});
        long startTime = System.nanoTime();
//        countHistogram格式不适用于该数据集集，换成构建pointCntList
        if (countHistogram.containsKey(i))
            countHistogram.put(i, countHistogram.get(i) + 1);
        else
            countHistogram.put(i, 1);
//        新的countHistogram，将过于离散的结果集中起来，也就是先增加一个集中的步骤
        int ii = i - i % 100;
        if (countHistogramNew.containsKey(ii))
            countHistogramNew.put(ii, countHistogramNew.get(ii) + 1);
        else
            countHistogramNew.put(ii, 1);
//        pointCntList.add(i);
        if (storeAllDatasetMemory) {
            System.out.println(fileNo);
            System.out.println(dataMapPorto.containsKey(fileNo));
            dataMapPorto.put(fileNo, a);
        }

        if (storeIndexMemory) {
            createDatasetIndex(fileNo, a);
        }
        if (!zcurveExist) {
            storeZcurve(a, fileNo);
        }
        long endTime = System.nanoTime();
        buildIndexTreeTime += (endTime - startTime) / 1000000000.0;
//		add createDatalake function so that EffectivenessStudy.main can run
//		createDatalake(1);
        //	return datasetMap;
    }

    public static void readContentBaidu(File file) throws IOException {
        System.out.println("read file " + file.getCanonicalPath());
//        long lineNumber = 0;
//        try (Stream<String> lines = Files.lines(file.toPath())) {
//            lineNumber = lines.count();
//        }
        List<double[]> l = new ArrayList<>();
//        double[][] a = new double[(int) lineNumber - 1][];// a little different
        int i = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String strLine;
            while ((strLine = br.readLine()) != null) {
                // System.out.println(strLine);
                String[] splitString = strLine.split(",");
                String aString = splitString[0];
                if (aString.equals("uid")) {
                    continue;
                }
                double[] b = new double[2];
                if (!splitString[12].equals("") && !splitString[11].equals("")) {
                    try {

//                        第一个值是纬度lat，第二个值是经度lng
                        b[0] = Double.parseDouble(splitString[12]);
                        b[1] = Double.parseDouble(splitString[11]);
                    } catch (Exception e) {
                        e.printStackTrace();
//                        System.out.println(e);
//                        System.out.println(file.getName());
//                        System.out.println(splitString[12]);
//                        System.out.println(splitString[11]);
                        continue;
                    }
//                    过滤掉不符合中国经纬度范围的数据
                    if (b[0] < 3 || b[0] > 54 || b[1] < 73 || b[1] > 136) {
                        continue;
                    }
                    l.add(b);
                    i++;
                }
            }
//			long endTime = System.nanoTime();
//			buildIndexTreeTime += (endTime - startTime) / 1000000000.0;
        }
        double[][] a = l.toArray(new double[][]{});
        long startTime = System.nanoTime();
//        countHistogram格式不适用于该数据集集，换成构建pointCntList
        if (countHistogram.containsKey(i))
            countHistogram.put(i, countHistogram.get(i) + 1);
        else
            countHistogram.put(i, 1);
//        新的countHistogram，将过于离散的结果集中起来，也就是先增加一个集中的步骤
        int ii = i - i % 100;
        if (countHistogramNew.containsKey(ii))
            countHistogramNew.put(ii, countHistogramNew.get(ii) + 1);
        else
            countHistogramNew.put(ii, 1);
//        pointCntList.add(i);
        if (storeAllDatasetMemory){
            System.out.println(fileNo);
            System.out.println(dataMapPorto.containsKey(fileNo));
            dataMapPorto.put(fileNo, a);
        }
        if (storeIndexMemory) {
            createDatasetIndex(fileNo, a);
        }
        if (!zcurveExist) {
            storeZcurve(a, fileNo);
        }
        long endTime = System.nanoTime();
        buildIndexTreeTime += (endTime - startTime) / 1000000000.0;
//		add createDatalake function so that EffectivenessStudy.main can run
//		createDatalake(1);
        //	return datasetMap;
    }


    /*
     * compute the normalization weight other dimensions based on the latitude (x) and longitude (y),
     * where the x and y indicate the position of latitude and longitude
     */
    static double[] normailizationWeight(Map<Integer, double[][]> datasetMap, int x, int y, int dimension, String filename) {
        double weight[] = new double[dimension];
        File tempFile = new File(filename);
        File tempFileRange = new File(filename + "range");
        if (tempFile.exists() && tempFileRange.exists()) {// if the serialization file exists, read the weights directly
            HashMap<Integer, Double> weightHashMap = EffectivenessStudy.deSerializationResult(filename);
            HashMap<Integer, Double> xyrange = EffectivenessStudy.deSerializationResult(filename + "range");
            for (int i = 0; i < dimension; i++) {
                weight[i] = weightHashMap.get(i);
            }
            minx = xyrange.get(1);
            miny = xyrange.get(2);
            spaceRange = xyrange.get(3);
        } else {
            HashMap<Integer, Double> weightHashMap = new HashMap<>();
            double max[] = new double[dimension];
            double min[] = new double[dimension];
            double range[] = new double[dimension];
            for (int i = 0; i < dimension; i++) {
                max[i] = Double.MIN_VALUE;
                min[i] = Double.MAX_VALUE;
            }
            for (int datasetid : datasetMap.keySet()) {
                double[][] dataset = datasetMap.get(datasetid);
                for (double[] datapoint : dataset) {
                    for (int i = 0; i < dimension; i++) {
                        if (datapoint[i] > max[i])
                            max[i] = datapoint[i];
                        if (datapoint[i] < min[i])
                            min[i] = datapoint[i];
                    }
                }
            }
            for (int i = 0; i < dimension; i++) {
                range[i] = max[i] - min[i];
                //	System.out.println(max[i]+","+min[i]+","+range[i]+","+datasetMap.keySet().size());
            }
            minx = min[x];
            miny = min[y];
            spaceRange = Math.max(range[x], range[y]);
            HashMap<Integer, Double> xyrange = new HashMap<>();
            xyrange.put(1, minx);
            xyrange.put(2, miny);
            xyrange.put(3, spaceRange);
            EffectivenessStudy.SerializedResultTable(filename + "range", xyrange);
            double scale = (range[x] + range[y]) / 2;
            for (int i = 0; i < dimension; i++) {
                if (i != x && i != y) {
                    weightHashMap.put(i, scale / range[i]);
                    weight[i] = scale / range[i];
                } else {
                    weightHashMap.put(i, 1.0);// storing the weights for normalization anywhere
                    weight[i] = 1;
                }
            }
            EffectivenessStudy.SerializedResultTable(filename, weightHashMap);    // no need to store another file but just the weight, serialize into mapping table.
        }
        return weight;// we will use the weight in the index building and search where distance computation is needed for radius
    }


    /*
     * read a folder and extract the corresponding column of each file inside
     */
    static Map<Integer, double[][]> readFolder(File folder, int limit) {
        File[] fileNames = folder.listFiles();
        int fileNo = 1;
        if (storeAllDatasetMemory)
            dataMapPorto = new HashMap<Integer, double[][]>();
        for (File file : fileNames) {
            if (file.isDirectory()) {
                readFolder(file, limit);
            } else {
                try {
                    String a = file.getName();
//            		if(a.length()<15){//argo
//            			datasetIdMapping.put(fileNo, a);
//            			readContent(file, fileNo++);
////            			add createDatalake function
////						createDatalake(1);
//            		}
//            		else {
//            			datasetIdMapping.put(fileNo, a);
//            			readContentChicago(file, fileNo++, a);
//            		}
//                    argoverse
                    if (a.length() < 15) {
                        datasetIdMapping.put(fileNo, a);
                        readContent(file, fileNo++);
                    } else {
                        datasetIdMapping.put(fileNo, a);
                        readContentChicago(file, fileNo++, a);
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            if (fileNo > limit)
                break;
        }
        return dataMapPorto;
    }

    /*
     * storing the z-curve
     */
    public static void storeZcurve(double[][] dataset, int datasetid) {
        if (zcodemap == null)
            zcodemap = new HashMap<Integer, ArrayList<Integer>>();
        int numberCells = (int) Math.pow(2, resolution);
        double unit = spaceRange / numberCells;
        ArrayList<Integer> zcodeaArrayList = new ArrayList<Integer>();
        for (int i = 0; i < dataset.length; i++) {
            int x = (int) ((dataset[i][0] - minx) / unit);
            int y = (int) ((dataset[i][1] - miny) / unit);
            int zcode = EffectivenessStudy.combine(x, y, resolution);
            if (!zcodeaArrayList.contains(zcode))
                zcodeaArrayList.add(zcode);
        }
        zcodemap.put(datasetid, zcodeaArrayList);
    }

    /*
     * read single dataset in an argoverse folder
     */
    public static void readContent(File file, int fileNo) throws IOException {
        //	System.out.println("read file " + file.getCanonicalPath());
        long lineNumber = 0;
        try (Stream<String> lines = Files.lines(file.toPath())) {
            lineNumber = lines.count();
        }
        double[][] a = new double[(int) lineNumber - 1][];// a little different
        int i = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String strLine;
            while ((strLine = br.readLine()) != null) {
                // System.out.println(strLine);
                String[] splitString = strLine.split(",");
//				3 for folder "data"
//				0 for folder "poi"
                String aString = splitString[3];
                if (aString.matches("-?\\d+(\\.\\d+)?")) {// only has float
                    a[i] = new double[2];
                    a[i][0] = Double.valueOf(splitString[3]);
                    a[i][1] = Double.valueOf(splitString[4]);
//					a[i][2] = Double.valueOf(splitString[0]);
//					a[i][0] = Double.valueOf(splitString[0]);
//					a[i][1] = Double.valueOf(splitString[1]);
                    i++;
                }
            }
//			long endTime = System.nanoTime();
//			buildIndexTreeTime += (endTime - startTime) / 1000000000.0;
        }
        long startTime = System.nanoTime();
        if (countHistogram.containsKey(i))
            countHistogram.put(i, countHistogram.get(i) + 1);
        else
            countHistogram.put(i, 1);

        //        新的countHistogram，将过于离散的结果集中起来，也就是先增加一个集中的步骤
        int ii = i - i % 10;
        if (countHistogramNew.containsKey(ii))
            countHistogramNew.put(ii, countHistogramNew.get(ii) + 1);
        else
            countHistogramNew.put(ii, 1);

        if (storeAllDatasetMemory)
            dataMapPorto.put(fileNo, a);
        if (storeIndexMemory) {
            createDatasetIndex(fileNo, a);
        }
        if (!zcurveExist) {
            storeZcurve(a, fileNo);
        }
        long endTime = System.nanoTime();
        buildIndexTreeTime += (endTime - startTime) / 1000000000.0;
//		add createDatalake function so that EffectivenessStudy.main can run
//		createDatalake(1);
        //	return datasetMap;
    }

    /*
     * read single dataset in a chicago folder
     */
    public static Map<Integer, double[][]> readContentChicago(File file, int fileNo, String filename) throws IOException {
        long lineNumber = 0;
        try (Stream<String> lines = Files.lines(file.toPath())) {
            lineNumber = lines.count();
        }
        System.out.println("file: " + filename);
//        double[][] a = new double[(int) lineNumber][];
        List<double[]> l = new ArrayList<>();
        int i = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String strLine;
            while ((strLine = br.readLine()) != null) {
                // System.out.println(strLine);
                String[] splitString = strLine.split(",");
				/*	Boolean isnumberBoolean = true;
				for(int j=0; j<11; j++)
					if(splitString[j].isEmpty() || splitString[j].matches("-?\\d+(\\.\\d+)?")==false) {
						isnumberBoolean = false;
						break;
					}
				if (isnumberBoolean) {*/
                //	write(filenameString, strLine+"\n");
                //	System.out.println(strLine);
//                b[i] = new double[dimension];
                try {
                    double[] b = new double[dimension];
                    b[0] = Double.parseDouble(splitString[7]);
                    b[1] = Double.parseDouble(splitString[8]);
                    b[2] = Double.parseDouble(splitString[9]);
                    b[3] = Double.parseDouble(splitString[10]);
                    b[4] = Double.parseDouble(splitString[0]);
                    b[5] = Double.parseDouble(splitString[1]);
                    b[6] = Double.parseDouble(splitString[2]);
                    b[7] = Double.parseDouble(splitString[3]);
                    b[8] = Double.parseDouble(splitString[4]);
                    b[9] = Double.parseDouble(splitString[5]);
                    b[10] = Double.parseDouble(splitString[6]);
                    l.add(b);
                    i++;
                } catch (NumberFormatException e) {
                    System.out.println("invalid number format");
//                    e.printStackTrace();
                }
            }
        }
        if (i > 0) {
            double[][] a = l.toArray(new double[][]{});
            if (countHistogram.containsKey(i))
                countHistogram.put(i, countHistogram.get(i) + 1);
            else
                countHistogram.put(i, 1);
            if (storeAllDatasetMemory)
                dataMapPorto.put(fileNo, a);
            if (storeIndexMemory) {
                createDatasetIndex(fileNo, a);
            }
            if (!zcurveExist)
                storeZcurve(a, fileNo);
        }
        return dataMapPorto;
    }

    /*
     * read from the mnist dataset and comparing with the dataset
     */
    public static Map<Integer, double[][]> readMnist(File file) throws IOException {
        Map<Integer, double[][]> datasetMap = new HashMap<Integer, double[][]>();
        //	System.out.println("read file " + file.getCanonicalPath());
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String strLine;
            int line = 1;
            while ((strLine = br.readLine()) != null) {
                // System.out.println(strLine);
                String[] splitString = strLine.split(" ");
                double[][] a = new double[28][];
                for (int i = 0; i < splitString.length; i++) {
                    if ((i + 1) % 28 == 0) {
                        a[i / 28] = new double[28];
                        for (int j = i - 27; j < i + 1; j++) {
                            //	System.out.println(j-i+27);
                            a[i / 28][j - i + 27] = Double.valueOf(splitString[j]);
                        }
                    }
                }
                datasetMap.put(line, a);
                line++;
            }
        }
        return datasetMap;
    }

    /*
     * read from the mnist dataset and comparing with the dataset
     */
    public static Map<Integer, double[][]> readMnist2D(File file) throws IOException {
        Map<Integer, double[][]> datasetMap = new HashMap<Integer, double[][]>();
        System.out.println("read file " + file.getCanonicalPath());
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String strLine;
            int line = 1;
            while ((strLine = br.readLine()) != null) {
                // System.out.println(strLine);
                String[] splitString = strLine.split(" ");
                double[][] a = new double[28][];
                for (int i = 0; i < splitString.length; i++) {
                    if ((i + 1) % 28 == 0) {
                        a[i / 28] = new double[28];
                        for (int j = i - 27; j < i + 1; j++) {
                            //	System.out.println(j-i+27);
                            a[i / 28][j - i + 27] = Double.valueOf(splitString[j]);
                        }
                    }
                }
                datasetMap.put(line, a);
                line++;
            }
        }
        return datasetMap;
    }

    /*
     * test whether hausdorff is metric
     */
    public void testTriangle() throws IOException {
        String foldername = "/Users/sw160/Desktop/argoverse-api/dataset/train/data";
        File folder = new File(foldername);
        Map<Integer, double[][]> dataMap = readFolder(folder, 1000);
        folder = new File("/Users/sw160/Downloads/torch-clus/dataset/minist_60000_784");
        dataMap = readMnist(folder);
        for (int a : dataMap.keySet()) {
            if (a > dataMap.size() - 3)
                break;
            //	System.out.println(dataMap.get(a)[0].length);
            double distance = Hausdorff.HausdorffExact(dataMap.get(1), dataMap.get(a), dataMap.get(a)[0].length, false);
            double distance1 = Hausdorff.HausdorffExact(dataMap.get(a), dataMap.get(a + 1), dataMap.get(a)[0].length, false);
            double distance2 = Hausdorff.HausdorffExact(dataMap.get(a + 1), dataMap.get(a + 2), dataMap.get(a)[0].length, false);
            double distance3 = Hausdorff.HausdorffExact(dataMap.get(a), dataMap.get(a + 2), dataMap.get(a)[0].length, false);
            if (distance1 + distance2 < distance3)
                System.out.println("no triangle inequality");//
        }
    }

    public static void testTriangleFrechet() throws IOException {
        String foldername = "/Users/sw160/Desktop/argoverse-api/dataset/train/data";
        File folder = new File(foldername);
        Map<Integer, double[][]> dataMap = readFolder(folder, 1000);
        folder = new File("/Users/sw160/Downloads/torch-clus/dataset/minist_60000_784");
        dataMap = readMnist(folder);
        for (int a : dataMap.keySet()) {
            if (a > dataMap.size() - 3)
                break;
            //	System.out.println(dataMap.get(a)[0].length);
            double distance = Hausdorff.Frechet(dataMap.get(1), dataMap.get(a), dataMap.get(a)[0].length);
            //	System.out.println(distance);
            double distance1 = Hausdorff.Frechet(dataMap.get(a), dataMap.get(a + 1), dataMap.get(a)[0].length);
            double distance2 = Hausdorff.Frechet(dataMap.get(a + 1), dataMap.get(a + 2), dataMap.get(a)[0].length);
            double distance3 = Hausdorff.Frechet(dataMap.get(a), dataMap.get(a + 2), dataMap.get(a)[0].length);
            if (distance1 + distance2 < distance3 || Math.abs(distance1 - distance2) > distance3)
                System.out.println("no triangle inequality");//
        }
    }

    public static void generateGTforPrediciton(String distanceString, Map<Integer, double[][]> dataMapPorto, Map<Integer, indexNode> indexMap, int dimension, int limit) {
        write(distanceString, "datasetID,dataset2ID,distance,normalizedDistance\n");
        for (int b : dataMapPorto.keySet()) {// we can build index to accelerate
            if (b > limit)
                break;
            for (int a : dataMapPorto.keySet()) {// we can build index to accelerate
                if (a > limit)
                    break;
                //		double distance = Hausdorff.HausdorffExact(dataMapPorto.get(b), dataMapPorto.get(a),
                //				dataMapPorto.get(a)[0].length, false);

                double distance = AdvancedHausdorff.IncrementalDistanceDirected(dataMapPorto.get(b), dataMapPorto.get(a), dimension, indexMap.get(b), indexMap.get(a), 0, 0, 0.05, null, null, null, true);
                // we should change it to our fast one, instead of this slow
                // building index for each dataset
                double pivot_distance = Util.EuclideanDis(indexMap.get(b).getPivot(), indexMap.get(a).getPivot(),
                        dimension);
                double max_distance = pivot_distance + indexMap.get(b).getRadius() + indexMap.get(a).getRadius();
                write(distanceString, Integer.toString(b) + "," + Integer.toString(a) + "," + Double.toString(distance) + ","
                        + Double.toString((max_distance - distance) / max_distance) + "\n");
            }
        }
    }

    public static void SerializedMappingTable(String file) {
        try {
            FileOutputStream fos = new FileOutputStream(file);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(datasetIdMapping);
            oos.close();
            fos.close();
            System.out.printf("Serialized HashMap data is saved in hashmap.ser");
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    /*
     * put the searched files into a folder
     */
    public static void UnionFile(String datasetid, String destination) throws FileNotFoundException, IOException {
        File file = new File(aString + "/" + datasetid);
        File dest = new File(destination + "/" + datasetid);
        try {
            FileUtils.copyFile(file, dest);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /*
     * this reads a single dataset after using the datalake index, instead of storing all the datasets in the main memory
     */
    public static double[][] readSingleFile(String datasetid) throws FileNotFoundException, IOException {
        File file = new File(aString + "/" + datasetid);
        long lineNumber = 0;
        try (Stream<String> lines = Files.lines(file.toPath())) {
            lineNumber = lines.count();
        }
        double[][] a = null;
        int i = 0;
        if (datalakeID == 4)
            try (BufferedReader br = new BufferedReader(new FileReader(file))) {// reading argo
                String strLine;
                a = new double[(int) lineNumber - 1][];
                while ((strLine = br.readLine()) != null) {
                    String[] splitString = strLine.split(",");
                    String aString = splitString[3];
                    if (aString.matches("-?\\d+(\\.\\d+)?")) {
                        a[i] = new double[3];
                        a[i][0] = Double.valueOf(splitString[3]);
                        a[i][1] = Double.valueOf(splitString[4]);
                        a[i][2] = Double.valueOf(splitString[0]);
                        i++;
                    }
                }
            }
        else if (datalakeID == 3)
            try (BufferedReader br = new BufferedReader(new FileReader(file))) {
                String strLine;
                a = new double[(int) lineNumber][];
                while ((strLine = br.readLine()) != null) {
                    String[] splitString = strLine.split(" ");
                    a[i] = new double[3];
                    a[i][0] = Double.valueOf(splitString[0]);
                    a[i][1] = Double.valueOf(splitString[1]);
                    a[i][2] = Double.valueOf(splitString[2]);
                    i++;
                }
            }
        else if (datalakeID == 5)
            try (BufferedReader br = new BufferedReader(new FileReader(file))) {
                String strLine;
                a = new double[(int) lineNumber][];
                while ((strLine = br.readLine()) != null) {
                    String[] splitString = strLine.split(",");
                    a[i] = new double[11];
                    a[i][0] = Double.valueOf(splitString[7]);
                    a[i][1] = Double.valueOf(splitString[8]);
                    a[i][2] = Double.valueOf(splitString[9]);
                    a[i][3] = Double.valueOf(splitString[10]);
                    a[i][4] = Double.valueOf(splitString[0]);
                    a[i][5] = Double.valueOf(splitString[1]);
                    a[i][6] = Double.valueOf(splitString[2]);
                    a[i][7] = Double.valueOf(splitString[3]);
                    a[i][8] = Double.valueOf(splitString[4]);
                    a[i][9] = Double.valueOf(splitString[5]);
                    a[i][10] = Double.valueOf(splitString[6]);
                    i++;
                }
            }
        return a;
    }

    /*
     * setting parameters and reading the whole datalake
     */
    static void readDatalake(int limit) throws IOException {
        File folder = new File(aString);
        countHistogram = new TreeMap<Integer, Integer>();
        //	dataMapPorto = new HashMap<>();
        long startTime = System.nanoTime();
//		valid option
        if (aString.contains("argo")) {// read datasets in a folder, such as argoverse and chicago
            datalakeID = 4;
//			dimension = 3;

            if (testFullDataset)
                limit = 270000;
            edgeString = "./index/dss/index/argo/argo_edge.txt";
            nodeString = "./index/dss/index/argo/argo_node.txt";
            zcodeSer = "./index/dss/index/argo/argo_bitmaps";
            distanceString = "./index/dss/index/argo/argo_haus_distance.txt";
            fetureString = "./index/dss/index/argo/argo_haus_features.txt";
            indexString = "./index/dss/index/argo/";
            dataMapPorto = readFolder(folder, limit + NumberQueryDataset);
//			add createDatalake function so that effectivenessStudy.main can run
//			createDatalake(1);
//			add SerializedMappingTable so that effectivenessStudy.main can run
//			SerializedMappingTable(indexString+"mapping.ser");
        }
//		valid option
        else if (aString.contains("chicago")) {// read datasets in a folder, such as argoverse and chicago
            datalakeID = 5;
            dimension = 11;
            if (testFullDataset)
//                limit = 280000;
                limit = 10000;
            edgeString = "./index/dss/index/chicago/chicago_edge.txt";
            nodeString = "./index/dss/index/chicago/chicago_node.txt";
            zcodeSer = "./index/dss/index/chicago/chicago_bitmaps";
            distanceString = "./index/dss/index/chicago/chicago_haus_distance.txt";
            fetureString = "./index/dss/index/chicago/chicago_haus_features.txt";
            indexString = "./index/dss/index/chicago/";
            dataMapPorto = readFolder(folder, limit + NumberQueryDataset);
            //	SerializedMappingTable(indexString+"mapping.ser");
        } else if (aString.contains("mnist")) {
            dimension = 28;
            edgeString = "./index/dss/mnist_edge.txt";
            nodeString = "./index/dss/mnist_node.txt";
            distanceString = "./index/dss/mnist_haus_distance.txt";
            fetureString = "./index/dss/mnist_haus_features.txt";
            dataMapPorto = readMnist(folder);
        } else if (aString.contains("beijing")) {
            dimension = 2;
            datalakeID = 2;
            if (testFullDataset)
                limit = 250000;
            edgeString = "./index/dss/index/beijing/beijing_edge.txt";
            nodeString = "./index/dss/index/beijing/beijing_node.txt";
            zcodeSer = "./index/dss/index/beijing/beijing_bitmaps";
            distanceString = "./index/dss/index/beijing/beijing_haus_distance.txt";
            fetureString = "./index/dss/index/beijing/beijing_haus_features.txt";
            indexString = "./index/dss/index/beijing/";
            dataMapPorto = readTrajectoryData(folder, limit);
        } else if (aString.contains("footprint")) {
            zcodeSer = "./index/footprint_bitmaps";
            edgeString = "./index/dss/nyc_edge.txt";
            nodeString = "./index/dss/nyc_node.txt";
            distanceString = "./index/dss/nyc_haus_distance.txt";
            fetureString = "./index/dss/nyc_haus_features.txt";
            dataMapPorto = readTrajectoryData(folder, limit);
//			valid option
        } else if (aString.contains("shapenet")) {
            datalakeID = 3;
            dimension = 3;
            if (testFullDataset) {
                limit = 10000;
            }
//                limit = 32000;
            zcodeSer = "./index/dss/index/shapenet/shapenet_bitmaps";
            edgeString = "./index/dss/index/shapenet/shapenet_edge.txt";
            nodeString = "./index/dss/index/shapenet/shapenet_node.txt";
            distanceString = "./index/dss/index/shapenet/shapenet_haus_distance.txt";
            fetureString = "./index/dss/index/shapenet/shapenet_haus_features.txt";
            indexString = "./index/dss/index/shapenet/";
            dataMapPorto = readShapeNet(folder, dataMapPorto, limit + NumberQueryDataset);
        } else if (aString.contains("porto")) {//reading porto
            datalakeID = 1;
            dimension = 2;
            if (testFullDataset)
                limit = 1000000;
            edgeString = "./index/dss/index/porto/porto_edge.txt";
            zcodeSer = "./index/dss/index/porto/porto_bitmaps";
            nodeString = "./index/dss/index/porto/porto_node.txt";
            distanceString = "./index/dss/index/porto/porto_haus_distance.txt";
            fetureString = "./index/dss/index/porto/porto_haus_features.txt";
            indexString = "./index/dss/index/porto/";
            dataMapPorto = readTrajectoryData(folder, limit);
        } else if (aString.contains("real")) {
            datalakeID = 6;
            dimension = 2;
            if (testFullDataset)
                limit = 3000;
            edgeString = "./index/dss/index/nyu/nyu_edge.txt";
            zcodeSer = "./index/dss/index/nyu/nyu_bitmaps";
            nodeString = "./index/dss/index/nyu/nyu_node.txt";
            distanceString = "./index/dss/index/nyu/nyu_haus_distance.txt";
            fetureString = "./index/dss/index/nyu/nyu_haus_features.txt";
            indexString = "./index/dss/index/nyu/";
            dataMapPorto = readNyu(folder, limit);
        } else if (aString.contains("btaa")) {
            datalakeID = 7;
            dimension = 2;
            if (testFullDataset)
                limit = 1000000;
            edgeString = "./index/dss/index/btaa/btaa_edge.txt";
            zcodeSer = "./index/dss/index/btaa/btaa_bitmaps";
            nodeString = "./index/dss/index/btaa/btaa_node.txt";
            distanceString = "./index/dss/index/btaa/btaa_haus_distance.txt";
            fetureString = "./index/dss/index/btaa/btaa_haus_features.txt";
            indexString = "./index/dss/index/btaa/";
            dataMapPorto = readNyu(folder, limit);
        }
        long endTime = System.nanoTime();
        double costTime = (endTime - startTime) / 1000000000.0;
        System.out.println("build index tree cost time: " + buildIndexTreeTime);
        System.out.println("readFolder cost time: " + costTime);
        System.out.println("dataset is loaded");
//		record covered points of each leaf indexNode
//		recordCoveredPointsNum();
//		record radius
        startTime = System.nanoTime();
//        recordRadius("radius_before");

//		startTime = System.nanoTime();
////		calculate max radius
//		calculateMaxRadius();
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("calculateMaxRadius cost time: " + costTime);

//		test
//		startTime = System.nanoTime();
//		double[][] data = calculateMaxRadiusByKneedleBefore();
//        calculateMaxRadiusByKneedle();
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("calculateMaxRadiusByKneedle cost time: " + costTime);
//
//		detect and delete outliers from index tree
//		startTime = System.nanoTime();
//		deleteOutlierFromAll();
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("deleteOutliers cost time: " + costTime);
//
//
//		update attributes
//		startTime = System.nanoTime();
//		updateIndexNodes();
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("updateIndexNodes cost time: " + costTime);
//
//        recordRadius("radius_after");
//		calculateMaxRadius();
//		calculateMaxRadiusByKneedleAfter(data);
//        calculateMaxRadiusByKneedle();
//        endTime = System.nanoTime();
//        deleteOutliersTime = (endTime - startTime) / 1000000000.0;
//        System.out.println("delete outliers cost time: " + deleteOutliersTime);
//        System.out.println("finish");

		/*startTime = System.nanoTime();
		createDatalake(1);
		endTime = System.nanoTime();
		costTime = (endTime - startTime) / 1000000000.0;
		System.out.println("createDatalake cost time: " + costTime);*/
        weight = normailizationWeight(dataMapPorto, 0, 1, dimension, indexString + "weight.ser");//the weight for each dimension
        clearBeforeWrite(indexString + "countHistogram.txt");
        for (int key : countHistogram.keySet()) {
            write(indexString + "countHistogram.txt", key + "," + countHistogram.get(key) + "\n");
        }
        clearBeforeWrite(indexString + "countHistogramConcentrate.txt");
        for (int key : countHistogramNew.keySet()) {
            write(indexString + "countHistogramConcentrate.txt", key + " " + countHistogramNew.get(key) + "\n");
        }
//        创建pointCntList，为了绘制数据分布直方图，专门为nyu和btaa数据集集设计
        if (pointCntList.size() > 0) {
            Collections.sort(pointCntList);
            int len = pointCntList.size();
            clearBeforeWrite(indexString + "countHistogramNew.txt");
            for (int i = 0; i < len; i++) {
                write(indexString + "countHistogramNew.txt", i + " " + pointCntList.get(i) + "\n");
            }
        }
    }

    //	record covered points of each leaf indexNode
    public static void recordCoveredPointsNum() {
        File file = new File(indexString + "coveredPointsNum.txt");
        if (datasetRoot != null) {
            if (coveredPointsMap.size() == 0) {
//				build map
                buildCoveredPointsMap(datasetRoot);
            }
            if (!file.exists()) {
                for (int key : coveredPointsMap.keySet()) {
                    write(indexString + "coveredPointsNum.txt", key + "," + coveredPointsMap.get(key) + "\n");
                }
            }
            System.out.println("Record covered points number!");
        } else {
            System.out.println("Cannot record covered points number!");
        }
    }

    //	build coveredPointsMap
    public static void buildCoveredPointsMap(indexNode node) {
        if (node.isLeaf()) {
            int coveredPoints = node.getTotalCoveredPoints();
//			test how covered points rise up to such large
            if (coveredPoints == 59) {
                System.out.println("?nb");
            }
            if (coveredPointsMap.containsKey(coveredPoints)) {
                coveredPointsMap.put(coveredPoints, coveredPointsMap.get(coveredPoints) + 1);
            } else {
                coveredPointsMap.put(coveredPoints, 1);
            }
        } else {
            for (indexNode childNode : node.getNodelist()) {
                buildCoveredPointsMap(childNode);
            }
        }
    }

    //	record radius of each leaf indexNode
    public static void recordRadius(String fileString) {
//		change objective from datasetRoot to indexNodes
        radiusList.clear();
        Stack<Double> radiusStack = new Stack<>();
        for (indexNode node : indexNodes) {
            buildRadiusQueue(node);
//			Collections.sort(radiusList);
        }
        while (!radiusQueue.isEmpty()) {
            double r = radiusQueue.poll();
            radiusList.add(r);
            radiusStack.push(r);
        }
        File file = new File(indexString + fileString + ".txt");
        if (file.exists()) {
            file.delete();
        }
        StringBuilder sb = new StringBuilder();
		/*for(double r: radiusList){
			sb.append(r).append("\n");
		}*/
        int size = radiusStack.size();
        for (int i = 0; i < size; i++) {
            sb.append(i).append(" ").append(radiusStack.pop()).append("\n");
        }
//		System.out.println(radiusStack.isEmpty());
        String str = new String(sb);
        String fileStr = indexString + fileString + ".txt";
        clearBeforeWrite(fileStr);
        write(fileStr, str);
        System.out.println("Record radius!");
    }

    //	build radiusSet
    public static void buildRadiusQueue(indexNode node) {
        if (node.isLeaf()) {
            double radius = node.getRadius();
//			radiusList.add(radius);
            radiusQueue.add(radius);
        } else {
            for (indexNode childNode : node.getNodelist()) {
                buildRadiusQueue(childNode);
            }
        }
    }

    public static boolean testRadiusAll(indexNode node, double[][] data) {
        boolean insideOfRange = true;
        insideOfRange = testRadius(node, data);
        if (!insideOfRange) {
            return false;
        }
        if (!node.isLeaf()) {
            for (indexNode child : node.getNodelist()) {
                insideOfRange = testRadius(child, data);
                if (!insideOfRange) {
                    break;
                }
                insideOfRange = testRadiusAll(child, data);
                if (!insideOfRange) {
                    break;
                }
            }
        }
        return insideOfRange;
    }

    public static boolean testRadius(indexNode node, double[][] data) {
        List<Integer> pointIdList = new ArrayList<>();
        node.getPointIdListAll(pointIdList);
        boolean insideOfRange = true;
        double[] pivot = node.getPivot();
        double radius = node.getRadius();
        double dist = 0;
        for (int i : pointIdList) {
            dist = Util.EuclideanDis(pivot, data[i], dimension);
            if (dist > radius) {
                System.out.println("distance from pivot to point " + i + " is " + dist +
                        ", larger than radius " + radius);
                insideOfRange = false;
                break;
            }
        }
        return insideOfRange;
    }

    //	record points in one leaf index node and write into file
    public static void recordPoints(double[][] dataset, List<Integer> pointIdList, double[] pivot) {
        if (dimension != 2) {
            return;
        }
        int size = pointIdList.size();
        double[][] points = new double[size][dimension];
        for (int i = 0; i < size; i++) {
//			has been changed
            points[i] = dataset[pointIdList.get(i)];
        }
        File file0 = new File(indexString + "points.txt");
        if (file0.exists()) {
            file0.delete();
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < size; i++) {
            sb.append(points[i][0]).append(", ").append(points[i][1]).append("\n");
        }
        String str0 = new String(sb);
        write(indexString + "points.txt", str0);

        File file1 = new File(indexString + "pivot.txt");
        if (file1.exists()) {
            file1.delete();
        }
        write(indexString + "pivot.txt", pivot[0] + ", " + pivot[1]);
        System.out.println("Record Points!");
    }

    //	calculate max radius by calling kneedle method
    public static double[][] calculateMaxRadiusByKneedleBefore() {
		/*Kneedle kn = new Kneedle();
		double[] x = new double[radiusList.size()];
		double[] y = new double[radiusList.size()];
		for(int i = 0; i < radiusList.size(); i++) {
			x[i] = i;
			y[i] = radiusList.get(i);
		}*/

        double[][] data = new double[radiusList.size()][2];
		/*SwingUtilities.invokeLater(() -> {
			PlotTest plot = new PlotTest(data);
			plot.setVisible(true);
			System.out.println("invoke");
		});*/
		/*PlotTest plot = new PlotTest(data);
		plot.setVisible(true);*/

        for (int i = 0; i < radiusList.size(); i++) {
            data[i][0] = i;
            data[i][1] = radiusList.get(i);
        }
        Kneedle kn = new Kneedle();
        ArrayList<double[]> result = kn.run(data, 30, 1, true);
//		System.out.println("there are " + result.size() + " knee points");
		/*for(double[] r: result) {
			for(double rr: r) {
				System.out.print(rr + " ");
			}
			System.out.println();
		}*/
        normalNodes = (int) result.get(result.size() - 1)[0];
        maxRadius = result.get(result.size() - 1)[1];
        System.out.println("new knee point:");
        System.out.println(normalNodes);
        System.out.println(maxRadius);

        return data;

//		int[] result = kn.detectElbowPoints(x, y);
		/*double rx = result[result.length - 1];
		double ry = radiusList.get((int)rx);
		System.out.println("[" + rx + ", " + ry + "]");*/
		/*for (int r: result) {
			System.out.println("[" + r + ", " + y[r] + "]");
		}*/
    }

    public static void calculateMaxRadiusByKneedleAfter(double[][] data) {
        double[][] dataNew = new double[radiusList.size()][2];

		/*PlotTest plot = new PlotTest(data);
		plot.setVisible(true);*/

        for (int i = 0; i < radiusList.size(); i++) {
            dataNew[i][0] = i;
            dataNew[i][1] = radiusList.get(i);
        }

		/*SwingUtilities.invokeLater(() -> {
			PlotTest plot = new PlotTest(dataNew, data);
			plot.setVisible(true);
		});*/

        Kneedle kn = new Kneedle();
        ArrayList<double[]> result = kn.run(dataNew, 30, 1, true);
//		System.out.println("there are " + result.size() + " knee points :");
		/*for(double[] r: result) {
			for(double rr: r) {
				System.out.print(rr + " ");
			}
			System.out.println();
		}*/
        normalNodes = (int) result.get(result.size() - 1)[0];
        maxRadius = result.get(result.size() - 1)[1];
        System.out.println("new new knee point:");
        System.out.println(normalNodes);
        System.out.println(maxRadius);
    }

    //	ultimate calculate max radius method
    public static void calculateMaxRadiusByKneedle() {
        int len = radiusList.size();
        double[][] data = new double[len][2];

        for (int i = 0; i < len; i++) {
            data[i][0] = i;
            data[i][1] = radiusList.get(i);
        }

        Kneedle kn = new Kneedle();
        ArrayList<double[]> result = kn.run(data, 30, 1, true);

        int length = result.size();
        normalNodes = (int) result.get(length - 1)[0];
        maxRadius = result.get(length - 1)[1];
        System.out.println("new knee point:");
        System.out.println(normalNodes);
        System.out.println(maxRadius);
    }

    //	find every node whose radius is too large and delete the outlier within
//	need to be divided into several functions
    public static void deleteOutlierFromAll() {
//		start from dataset rootNode rather than datalake rootNode
        int len = indexNodes.size();
        for (int i = 0; i < len; i++) {
            double[][] dataset = dataMapPorto.get(i + 1);
            indexNode node = indexNodes.get(i);
            deleteOutlierFromDataset(dataset, node, i);
        }
    }

    //	detect and delete outliers of one dataset
    public static void deleteOutlierFromDataset(double[][] dataset, indexNode node, int i) {
        if (node.isLeaf() && node.getRadius() > maxRadius) {
//			for test
			/*System.out.println("dataset " + i + " has outliers");
			System.out.println("original totalCoveredPoints is " + node.getTotalCoveredPoints());
			System.out.println("original radius is " + node.getRadius());*/
            deleteOutlierFromNode(dataset, node);
//			for test
			/*System.out.println("original totalCoveredPoints is " + node.getTotalCoveredPoints());
			System.out.println("original radius is " + node.getRadius());*/
        } else if (!node.isLeaf()) {
            for (indexNode child : node.getNodelist()) {
                deleteOutlierFromDataset(dataset, child, i);
            }
        }
    }

    //	detect and delete outliers of one indexNode
    public static void deleteOutlierFromNode(double[][] dataset, indexNode node) {
//		recalculate pivot
//		node.setPivot(dataset);
//		set new pivot
//		double[] pivot = node.getPivot();
        int coveredPoints = node.getTotalCoveredPoints();
//		get distance of new pivot to each point
//		if larger than max radius, delete the point
        node.sortPoints(dataset);
        Iterator<Integer> it = node.getpointIdList().iterator();
//		boolean firstSet = true;
        boolean isChanged = false;
        int before = node.getpointIdList().size();
        if (pointsRecorded == -1 && before == 19) {
            recordPoints(dataset, node.getpointIdList(), node.getPivot());
            pointsRecorded++;
        }
        while (it.hasNext()) {
//			double[] point = dataset[it.next() - 1];
            double[] point = dataset[it.next()];
//			call an Euclidean distance calculation method
//			double radiusNew = Util.EuclideanDis(pivotNew, point, dimension);
            double radiusNew = node.getPivotdistance(point);
            if (radiusNew > maxRadius) {
                it.remove();
                coveredPoints--;
                node.setPivotAndSum(dataset);
                isChanged = true;
            }
			/*if (firstSet) {
				node.setPivotAndSum(dataset);
				firstSet = false;
			}*/
        }
//		if(isChanged) {
//			System.out.println(before + " - " + node.getpointIdList().size() + " = "
//			+ (before - node.getpointIdList().size()));
//		}
//		recalculate pivot and radius
//		node.setPivotAndSum(dataset);
        node.setRadius(dataset);

//		renew other attributes
        node.setTotalCoveredPoints(coveredPoints);

        if (pointsRecorded == 0) {
            recordPoints(dataset, node.getpointIdList(), node.getPivot());
            pointsRecorded++;
        }
    }


    //	update all indexNodes' attributes
    public static void updateIndexNodes() {
//		traverse indexNodes and go into every dataset rootNode
//		take a postorder traversal into every dataset rootNode
//		update pivot ?, radius ?, totalCoveredPoints ?, sum ?, mbrmax ?, mbrmin ?, signature
        int len = indexNodes.size();
        for (int i = 0; i < len; i++) {
            double[][] dataset = dataMapPorto.get(i + 1);
            indexNode node = indexNodes.get(i);
            updateIndexNode(dataset, node);
        }
    }

    //	update one dataset indexNode's attributes
    public static void updateIndexNode(double[][] dataset, indexNode node) {
        indexAlgorithm<indexNode> ia = new indexAlgorithm<>();
        if (!node.isLeaf()) {
            for (indexNode childNode : node.getNodelist()) {
                updateIndexNode(dataset, childNode);
            }
            double[] pivot = new double[dimension];
            double[] sum = new double[dimension];
            double radius = 0;
            int totalCoveredPoints = 0;
            for (indexNode child : node.getNodelist()) {
                for (int i = 0; i < dimension; i++) {
                    sum[i] += child.getSum()[i];
                }
                totalCoveredPoints += child.getTotalCoveredPoints();
            }
            for (int i = 0; i < dimension; i++) {
//				pivot[i] = sum[i] / node.getNodelist().size();
                pivot[i] = sum[i] / totalCoveredPoints;
            }
            for (indexNode child : node.getNodelist()) {
                radius = Math.max(radius, Util.EuclideanDis(pivot, child.getPivot(), dimension) + child.getRadius());
            }
            ia.calculateMaxBoundBox(node, dimension, dataset);
            ia.calculateMinBoundBox(node, dimension, dataset);
            node.setPivot(pivot);
            node.setRadius(radius);
            node.setSum(sum);
            node.setTotalCoveredPoints(totalCoveredPoints);
        }
    }

    //	call python program to calculate the knee point as the max radius
    public static void calculateMaxRadius() {
        Process proc;
        try {
//			proc = Runtime.getRuntime().exec("python \"D:\\Projects\\Spadas\\kneed_test.py\"");
            proc = Runtime.getRuntime().exec("python \".\\src\\kneed_test.py\"");
            BufferedReader in = new BufferedReader(new InputStreamReader(proc.getInputStream()));
//			System.out.println(in.readLine());
//			System.out.println(in.readLine());
            normalNodes = Integer.parseInt(in.readLine());
            maxRadius = Double.parseDouble(in.readLine());
            System.out.println(normalNodes);
            System.out.println(maxRadius);
            in.close();
            proc.waitFor();
        } catch (InterruptedException | IOException e) {
            e.printStackTrace();
        }
    }

    // write the information into files
    public static void write(String fileName, String content) {
        RandomAccessFile randomFile = null;
//        System.out.println("content =");
//        System.out.println(content);
        try {
            randomFile = new RandomAccessFile(fileName, "rw");
            long fileLength = randomFile.length();
            randomFile.seek(fileLength);
//            randomFile.setLength(0);
//            randomFile.seek(0);
            randomFile.writeBytes(content);
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (randomFile != null) {
                try {
                    randomFile.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public static void clearBeforeWrite(String fileName) {
        RandomAccessFile randomFile = null;
        try {
            randomFile = new RandomAccessFile(fileName, "rw");
//            long fileLength = randomFile.length();
//            randomFile.seek(0);
            randomFile.setLength(0);
//            randomFile.seek(0);
//            randomFile.writeBytes(content);
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (randomFile != null) {
                try {
                    randomFile.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    // create index for a single dataset in a datalake
    static void createDatasetIndex(int a, double[][] dataset) {
        indexNode rootBall;
        if (buildOnlyRoots) {// just create the root node, for datalake creation
            rootBall = createRootsDataset(dataset, dimension, a);
        } else {// create the full tree
            long startTime = System.nanoTime();
            rootBall = indexDSS.buildBalltree2(dataset, dimension, capacity, null, null, weight);
            long endTime = System.nanoTime();
            double costTime = (endTime - startTime) / 1000000000.0;
//			System.out.println("buildBalltree cost time: " + costTime);
            String indexFileName = indexString + String.valueOf(a) + ".txt";
            File tempFile = new File(indexFileName);
            if (!tempFile.exists() && saveDatasetIndex)//&& !aString.contains("porto") && !aString.contains("beijing")
                indexDSS.storeIndex(rootBall, 1, indexFileName, 0);
            if (indexMap == null)
                indexMap = new HashMap<>();
            indexMap.put(a, rootBall);
        }
        indexDSS.setGloabalid();
        rootBall.setroot(a);// set an id to identify which dataset it belongs to
        if (a < limit)
            indexNodes.add(rootBall);
    }

    /*
     * just create a root node to cover all the points
     */
    static indexNode createRootsDataset(double[][] dataset, int dimension, int datasetid) {
        double max[], min[], pivot[], radius = 0;
        max = new double[dimension];
        min = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            max[i] = Double.MIN_VALUE;
            min[i] = Double.MAX_VALUE;
        }
        pivot = new double[dimension];
        for (double[] datapoint : dataset) {
            for (int i = 0; i < dimension; i++) {
                if (datapoint[i] > max[i])
                    max[i] = datapoint[i];
                if (datapoint[i] < min[i])
                    min[i] = datapoint[i];
            }
        }
        for (int i = 0; i < dimension; i++) {
            pivot[i] = (max[i] + min[i]) / 2;
            if (weight != null)
                radius += Math.pow((max[i] - min[i]) * weight[i] / 2, 2); // the radius should consider d dimension
            else
                radius += Math.pow((max[i] - min[i]) / 2, 2);
        }
        radius = Math.sqrt(radius);
        indexNode rootNode = new indexNode(dimension);
        rootNode.setMBRmax(max);
        rootNode.setMBRmax(min);
        rootNode.setRadius(radius);
        rootNode.setPivot(pivot);
        rootNode.setTotalCoveredPoints(dataset.length);
        return rootNode;
    }

    /*
     * get the count of nodes in the lightweight tree,
     */
    static long getDatalakeIndexSize(indexNode root, boolean centorid) {
        if (root == null)
            return 0;
        if (root.isrootLeaf()) {
            return root.getMemorySpadas(dimension, centorid);
        } else {
            Set<indexNode> listnode = root.getNodelist();
            long max = 0;
            for (indexNode aIndexNode : listnode) {
                max += aIndexNode.getMemorySpadas(dimension, centorid);
                max += getDatalakeIndexSize(aIndexNode, centorid);
            }
            return max;
        }
    }

    /*
     * get the count of nodes in the lightweight tree,
     */
    static long getDatasetIndexSize(indexNode root, boolean centorid) {
        if (root == null)
            return 0;
        if (root.isLeaf()) {
            return root.getMemorySpadas(dimension, centorid);
        } else {
            Set<indexNode> listnode = root.getNodelist();
            long max = 0;
            for (indexNode aIndexNode : listnode) {
                max += aIndexNode.getMemorySpadas(dimension, centorid);
                max += getDatasetIndexSize(aIndexNode, centorid);
            }
            return max;
        }
    }

    /*
     * get the count of nodes in the lightweight tree,
     */
    static long getDatalakeIndexSizeDisk(int dimension) {
        if (datalakeIndex == null)
            return 0;
        else {
            return datalakeIndex.size() * datalakeIndex.get(1).getMemorySpadas(dimension, null);
        }
    }

    static long getDatalakeSize(int dimension) {
        if (dataMapPorto == null) {
            return 0;
        } else {
            long size = 0;
            for (int datasetid : dataMapPorto.keySet()) {
                double[][] dataset = dataMapPorto.get(datasetid);
                size += dataset.length * dimension * 8;
            }
            return size;
        }
    }

    /*
     * get the count of nodes in the lightweight tree,
     */
    static long getDatasetIndexSizeDisk() {
        if (datasetIndex == null)
            return 0;
        else {
            long size = 0;
            for (int datasetid : datasetIndex.keySet()) {
                Map<Integer, indexNode> nodeMap = datasetIndex.get(datasetid);
                size += nodeMap.size() * nodeMap.get(1).getMemorySpadas(dimension, null);
            }
            return size;
        }
    }

    static long getAllDatasetIndexSize() {
        long max = 0;
        if (indexMap != null)
            for (int i : indexMap.keySet()) {
                indexNode root = indexMap.get(i);
                max += getDatasetIndexSize(root, false);
            }
        return max;
    }

    static void serializeIndexNode(indexNode node) {
        System.out.println(node.getDatasetID());
        String filePath = indexString + "index_node.dat";
        File file = new File(filePath);
        if (!file.exists()) {
            try {
                file.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        try (ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(filePath))) {
            oos.writeObject(node);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static indexNode deserializeIndexNode(String filePath) {
        try (ObjectInputStream ois = new ObjectInputStream(new FileInputStream(filePath))) {
            datasetRoot = (indexNode) ois.readObject();
        } catch (IOException | ClassNotFoundException e) {
            e.printStackTrace();
        }
        System.out.println(datasetRoot.getDatasetID());
        return datasetRoot;
    }


    /*
     * generate the features for index
     */
    static void generateIndexFeatures() throws IOException {
        double[][] dataset = new double[1917494][300];
        edgeString = "./index/dss/glove_edge.txt";
        nodeString = "./index/dss/glove_node.txt";
        String fileString = "/Users/sw160/Desktop/glove.42B.300d.txt";
        File file = new File(fileString);
        int count = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String strLine;
            while ((strLine = br.readLine()) != null) {
                String[] splitString = strLine.split(" ");
                //	System.out.println(strLine);
                //	System.out.println(splitString.length);
                if (splitString.length <= 300)
                    continue;
                for (int i = 1; i <= 300; i++) {
                    dataset[count][i - 1] = Double.valueOf(splitString[i]);
                }
                count++;
            }
        }
        dimension = 300;
        capacity = 30;
        indexNode rootBall = indexDSS.buildBalltree2(dataset, dimension, capacity, null, null, weight);
        write(edgeString, "datasetID,startNode,endNode,distance\n");
        write(nodeString, "datasetID,node,#CoverPoints,Depth,Radius,Leaf,PivotPoint\n");
        indexDSS.storeFeatures(rootBall, 1, 1, edgeString, nodeString, 1);
        // get the features of index tree,
        //	double features[] = datasetFeatures.getFeature(fetureString, rootBall,rootBall.getTotalCoveredPoints(), capacity, a);
        // featureMap.put(a, features);
    }

    // create datalake index based on the root nodes of all datasets, too slow
//    建立上层索引
    static void createDatalake(int N) {
        indexDSS.setGloabalid();
//        如果下层索引已经在内存，已经在indexNodes数据结构中
//        如果需要将上层索引序列化到disk中，则重新创建上层索引并序列化到disk中
        if (storeIndexMemory) {// build index for the whole datalake
            datasetRoot = indexDSS.indexDatasetKD(indexNodes, dimension, capacity, weight);
            System.out.println("index created");
            //	if(storeIndex)
//            先清除文件中内容
            clearBeforeWrite(indexString + "datalake" + N + "-" + capacity + "-" + dimension + ".txt");
//            再将上层索引写入文件
//            不是这索引有啥用啊
//            太耗时间了先注释掉，需要用再注释回来
//            indexDSS.storeDatalakeIndex(datasetRoot, 1, indexString + "datalake" + N + "-" + capacity + "-" + dimension + ".txt", 0);//store the
//            试图序列化全部索引，但耗时过长，占用内存过大，只能放弃
//            serializeIndexNode(datasetRoot);
            //	datalakeIndex = new HashMap<>();
            //	indexDSS.reorgnizeIndex(datasetRoot, 1, 0, datalakeIndex);//putting all nodes into hashmap
        } else {// 否则下层索引不在内存，不存在indexNodes数据结构来保存节点
            datalakeIndex = indexDSS.restoreDatalakeIndex(indexString + "datalake" + N + "-" + capacity + "-" + dimension + ".txt", dimension);//the datalake index
//            搞完只有上层索引没有下层索引了，不知道能干嘛
            datasetRoot = deserializeIndexNode(indexString + "index_node.dat");
        }

        if (datalakeIndex != null) {
            datalakeIndex.get(1).setMaxCovertpoint(datalakeIndex); // for the query that measures how many points in the intersected range
            datalakeIndex.get(1).buildsignarture(zcodemap, datalakeIndex); // for the gird-based overlap query, need to debug for the memory mode
        } else {
            datasetRoot.setMaxCovertpoint(datalakeIndex);
            datasetRoot.buildsignarture(zcodemap, datalakeIndex);
        }
    }

    // prepare two large scale datasets for pairwise comparison
    static Map<Integer, double[][]> prepareLargeDatasets(int numberofDatasets) {
        Map<Integer, double[][]> dataMapLargeMap = new HashMap<>();
        int count = 0;
        int length = 0, length1 = 0, length2 = 0;
        for (int id : dataMapPorto.keySet()) {
            double[][] dataset = dataMapPorto.get(id);
            length += dataset.length;
            count++;
            if (count == numberofDatasets) {
                length1 = length;
                length = 0;
            } else if (count >= 2 * numberofDatasets) {
                length2 = length;
                break;
            }
        }
        double[][] newdataset = new double[length1][dimension];
        double[][] newdataset1 = new double[length2][dimension];
        count = 0;
        length1 = 0;
        length2 = 0;
        for (int id : dataMapPorto.keySet()) {
            double[][] dataset = dataMapPorto.get(id);
            count++;
            if (count <= numberofDatasets) {
                for (int i = 0; i < dataset.length; i++) {
                    for (int j = 0; j < dimension; j++) {
                        newdataset[i + length1][j] = dataset[i][j];
                    }
                }
                length1 += dataset.length;
            } else if (count <= 2 * numberofDatasets && count > numberofDatasets) {
                for (int i = 0; i < dataset.length; i++) {
                    for (int j = 0; j < dimension; j++) {
                        newdataset1[i + length2][j] = dataset[i][j];
                    }
                }
                length2 += dataset.length;
            } else {
                break;
            }
        }
        dataMapLargeMap.put(1, newdataset);
        dataMapLargeMap.put(2, newdataset1);
        return dataMapLargeMap;
    }

    /*
     * four methods, existing methods: R-tree, Pami, ball-tree, half-ball theory
     */
    static Pair<Double, PriorityQueue<queueMain>> testPairwiseHausdorff(int queryID, int datasetID, int testOrder, int timeID, boolean fastOnly) throws FileNotFoundException, IOException {
        // exact computation
        numberofComparisons = 4;
        int countAlgorithm = 0;
        double distance = 0;

        System.out.println("*******testing pairwise Hausdorff");
        // brute force, too slow, we ignore it here
        long startTime1 = System.nanoTime();
//		distance = Hausdorff.HausdorffExact(dataMapPorto.get(queryID), dataMapPorto.get(datasetID), dataMapPorto.get(datasetID)[0].length, false);
        System.out.println("exact hausdorff distance: " + distance);
        long endtime = System.nanoTime();
        //	System.out.println((endtime-startTime1)/1000000000.0);
        //	System.out.println(distance);
        //	timeMultiple[0][testOrder*numberofComparisons+countAlgorithm++] = (endtime-startTime1)/1000000000.0;
        indexNode queryNode = null, datanode = null;
        Map<Integer, indexNode> queryindexmap = null, dataindexMap = null; // datasetIndex.get(queryid);
        double[][] querydata, dataset;
        if (dataMapPorto == null) {
            querydata = readSingleFile(datasetIdMapping.get(queryID));
            dataset = readSingleFile(datasetIdMapping.get(datasetID));
        } else {
            querydata = dataMapPorto.get(queryID);
            dataset = dataMapPorto.get(datasetID);
        }
        if (indexMap != null) {
            queryNode = indexMap.get(queryID);
            datanode = indexMap.get(datasetID);
        } else {
            if (datasetIndex != null) {
                queryindexmap = datasetIndex.get(queryID);
                dataindexMap = datasetIndex.get(datasetID);
                queryNode = queryindexmap.get(1);
                datanode = dataindexMap.get(1);
            } else {
                queryNode = indexDSS.buildBalltree2(querydata, dimension, capacity, null, null, weight);
                datanode = indexDSS.buildBalltree2(dataset, dimension, capacity, null, null, weight);
            }
        }

        // 1, using knn to search
        startTime1 = System.nanoTime();
        if (!fastOnly) {
            double pivot_distance = Util.EuclideanDis(queryNode.getPivot(), datanode.getPivot(), dimension);
            double ub = pivot_distance + queryNode.getRadius() + datanode.getRadius();
            distance = Hausdorff.HausdorffWithIndexDFS(querydata, dataset, dimension,
                    queryNode, datanode, 0, ub, indexDSS);
        }
        endtime = System.nanoTime();
        //	System.out.println(distance);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;


        // 2, pvldb 11 using MBR bounds
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        AdvancedHausdorff.setBoundChoice(2); // setting the bound option
        if (!fastOnly) {
            AdvancedHausdorff.IncrementalDistance(querydata, dataset, dimension, queryNode, datanode, 0, 1,
                    0, false, 0, false, queryindexmap, dataindexMap, null,
                    true);
        }
        endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0 + "," + AdvancedHausdorff.disCompTime + ","
                + AdvancedHausdorff.EuclideanCounter);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;

        // 3, pami 15
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        if (!fastOnly) {
            if (querydata.length < 10000)// will be out of memory if not
                distance = Hausdorff.earlyBreaking(querydata, dataset, dimension, false);
        }
        endtime = System.nanoTime();

        System.out.println((endtime - startTime1) / 1000000000.0);
        //	System.out.println(distance);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;

        // 4, our fast algorithm
        AdvancedHausdorff.setBoundChoice(0);
        startTime1 = System.nanoTime();//exact search: 0, 1, 0
        Pair<Double, PriorityQueue<queueMain>> resultPair = AdvancedHausdorff.IncrementalDistance(querydata, dataset, dimension, queryNode, datanode, 0, 1, 0, false, 0, false, queryindexmap, dataindexMap, null, true);
        endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0 + "," + AdvancedHausdorff.disCompTime + "," + AdvancedHausdorff.EuclideanCounter);
        //	System.out.println(resultPair.getLeft());
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;
        return resultPair;

		/*
		 *
		// the basic version of our method
		startTime1 = System.nanoTime();
	//	distance = Hausdorff.IncrementalDistance(dataMapPorto.get(queryID), dataMapPorto.get(a), dimension, indexMap.get(1), indexMap.get(a), 1);
		distance = AdvancedHausdorff.IncrementalDistance(dataMapPorto.get(queryID), dataMapPorto.get(datasetID), dimension, indexMap.get(queryID), indexMap.get(datasetID), 0);
		endtime = System.nanoTime();
		System.out.println((endtime-startTime1)/1000000000.0);
		System.out.println(distance);
		timeMultiple[0][testOrder*numberofComparisons+countAlgorithm++] = (endtime-startTime1)/1000000000.0;


		startTime1 = System.nanoTime();//exact search: 0, 1, 0
		distance = AdvancedHausdorff.IncrementalDistanceDirected(dataMapPorto.get(queryID), dataMapPorto.get(datasetID), dimension, indexMap.get(queryID), indexMap.get(datasetID), 0, 1, 0, null, null, null, true);
		endtime = System.nanoTime();
		System.out.println((endtime-startTime1)/1000000000.0+ ","+ AdvancedHausdorff.disCompTime+","+AdvancedHausdorff.EuclideanCounter);
		System.out.println(distance);//+","+resultPair1.getLeft()+","+resultPair2.getLeft()
		System.out.println();
		timeMultiple[0][testOrder*numberofComparisons+countAlgorithm++] = (endtime-startTime1)/1000000000.0;
		*/
    }

    /*
     * we test the range query for union operation, just pairise
     */
    static void testPairwiseRangeQuery(indexNode dataseNode, double[][] dataset, int testOrder, int timeID) {
        numberofComparisons = 4;

        int countAlgorithm = 3;

        double[] querymax = new double[dimension];
        double[] querymin = new double[dimension];
        for (int i = 0; i < dimension; i++) {// to enlarge the range for evaluation
            querymax[i] = dataseNode.getPivot()[i] + error;
            querymin[i] = dataseNode.getPivot()[i] - error;
        }

        // baseline method here, scan every datasets

        long startTime1 = System.nanoTime();
        double min_dis = Double.MAX_VALUE;
        ArrayList<double[]> points = null;
        dataseNode.coveredPOints(querymax, querymin, min_dis, dimension, dataset, null, dimNonSelected, dimensionAll, points);
        long endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;
    }

    // restored index, not putting all the datasets and index in memeory, only access when it is checked after passing datalake index
    static void testRestoredIndex() {
        /* test the restored index, it works well, while it is a little slower because of finding the hashmap */
        long startTime1 = System.nanoTime();//exact search: 0, 1, 0
        Pair<Double, PriorityQueue<queueMain>> resultPairq = AdvancedHausdorff.IncrementalDistance(dataMapPorto.get(1), dataMapPorto.get(3), dimension, indexMap.get(1), indexMap.get(3), 0, 1, 0, false, 0, false, null, null, null, true);
        long endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0 + "," + AdvancedHausdorff.disCompTime + "," + AdvancedHausdorff.EuclideanCounter);
        System.out.println(resultPairq.getLeft());

        startTime1 = System.nanoTime();//exact search: 0, 1, 0
        datasetIndex = indexDSS.restoreIndex(indexString, dimension, dataMapPorto);
        endtime = System.nanoTime();
        //	System.out.println("building index costs: "+(endtimea-startTime1a)/1000000000.0);
        System.out.println("reconstruction costs: " + (endtime - startTime1) / 1000000000.0);

        startTime1 = System.nanoTime();//exact search: 0, 1, 0
        Pair<Double, PriorityQueue<queueMain>> resultPair11 = AdvancedHausdorff.IncrementalDistance(dataMapPorto.get(1), dataMapPorto.get(3), dimension, datasetIndex.get(1).get(1), datasetIndex.get(3).get(1), 0, 1, 0, false, 0, false, datasetIndex.get(1), datasetIndex.get(3), null, true);
        endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0 + "," + AdvancedHausdorff.disCompTime + "," + AdvancedHausdorff.EuclideanCounter);
        System.out.println(resultPair11.getLeft());
    }

    static indexNode getQueryNode(double[][] querydata, int queryid) throws FileNotFoundException, IOException {
        indexNode queryNode = null;
        if (dataMapPorto == null)
            querydata = readSingleFile(datasetIdMapping.get(queryid));
        else
            querydata = dataMapPorto.get(queryid);
        if (indexMap != null) {
            System.out.println();
            queryNode = indexMap.get(queryid);
        } else {
            if (datasetIndex != null) {
                queryNode = datasetIndex.get(queryid).get(1);
            } else {
                queryNode = indexDSS.buildBalltree2(querydata, dimension, capacity, null, null, weight);
            }
        }
        return queryNode;
    }

    // four methods: scanning-all, sigspatial10, early adandoning, datalake
    static void testTopkHausdorff(int k, int testOrder, int queryid, int timeID, int limit) throws FileNotFoundException, IOException {
        numberofComparisons = 4;
        int countAlgorithm = 0;
        /* read query dataset */
        indexNode queryNode = null;
        Map<Integer, indexNode> queryindexmap = null; // datasetIndex.get(queryid);
        double[][] querydata;
        if (dataMapPorto == null)
            querydata = readSingleFile(datasetIdMapping.get(queryid));
        else
            querydata = dataMapPorto.get(queryid);
        if (indexMap != null) {
            queryNode = indexMap.get(queryid);
        } else {
            if (datasetIndex != null) {
                queryindexmap = datasetIndex.get(queryid);
                queryNode = queryindexmap.get(1);
            } else {
                queryNode = indexDSS.buildBalltree2(querydata, dimension, capacity, null, null, weight);
            }
        }

        System.out.println("*******testing top-k Hausdorff");
        AdvancedHausdorff.setBoundChoice(0);//using our fast ball bound

        // 1, conduct top-k dataset search, the baseline, scan one by one
        long startTime1 = System.nanoTime();// exact search: 0, 1, 0
        		Search.hausdorffEarylyAbandonging(dataMapPorto, queryid, indexMap, dimension, limit, false, saveDatasetIndex,
                        dimNonSelected, dimensionAll, querydata, datasetIdMapping, capacity, weight, queryindexmap, queryNode, indexString);//no top-k early breaking
        long endtime = System.nanoTime();
        System.out.println("top-k search costs: " + (endtime - startTime1) / 1000000000.0);
//        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;


        // 2, early abandoning, compute the bound first, then see whether to abandon it
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        Search.hausdorffEarylyAbandonging(dataMapPorto, queryid, indexMap, dimension, limit, true, saveDatasetIndex,
                dimNonSelected, dimensionAll, querydata, datasetIdMapping, capacity, weight, queryindexmap, queryNode, indexString);
        endtime = System.nanoTime();
        System.out.println("top-k search costs: " + (endtime - startTime1) / 1000000000.0);
//        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;

        // 3, rank all the candidate by lower bounds, then one by one
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        int datasetID = Search.HausdorffEarylyAbandongingRanking(dataMapPorto, queryid, indexMap, dimension, limit,
                resultPair, queryindexmap, saveDatasetIndex, dimNonSelected, dimensionAll, error, querydata,
                datasetIdMapping, queryNode, capacity, weight, indexString);
        endtime = System.nanoTime();
        System.out.println("top-k search costs: " + (endtime - startTime1) / 1000000000.0);
//        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;

        // 4, using datalake index to prune
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        HashMap<Integer, Double> result = Search.pruneByIndex(dataMapPorto, datasetRoot, queryNode, queryid,
                dimension, indexMap, datasetIndex, queryindexmap, datalakeIndex, datasetIdMapping, k,
                indexString, dimNonSelected, dimensionAll, error, capacity, weight, saveDatasetIndex, querydata);
        endtime = System.nanoTime();
        System.out.println("top-k search costs: " + (endtime - startTime1) / 1000000000.0);
        System.out.println();
//        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;
    }

    /*
     * test the approximate version of knn and pairwise by increasing the error,
     * we have two comparisons, pairwise and Hausdorff
     */
    static void testHausdorffApproxi(int k, int testOrder, int queryid, int timeID) throws FileNotFoundException, IOException {
        numberofComparisons = 4;
        int countAlgorithm = 2;

        System.out.println("*******testing approximate Hausdorff");

        // 1, pairwise
        int datasetID = queryid + 1;
        long startTime1 = System.nanoTime();
        double dataMap[][], query[][];
        if (dataMapPorto != null) {
            dataMap = dataMapPorto.get(datasetID);
            query = dataMapPorto.get(queryid);
        } else {
            dataMap = Framework.readSingleFile(datasetIdMapping.get(datasetID));
            query = Framework.readSingleFile(datasetIdMapping.get(queryid));
        }
        indexNode queryNode, datanode;
        Map<Integer, indexNode> queryindexmap = null, dataindexMap = null;
        if (indexMap != null) {
            queryNode = indexMap.get(queryid);
            datanode = indexMap.get(datasetID);
        } else {
            if (datasetIndex != null) {
                queryindexmap = datasetIndex.get(queryid);
                dataindexMap = datasetIndex.get(datasetID);
                queryNode = queryindexmap.get(1);
                datanode = dataindexMap.get(1);
            } else {
                queryNode = indexDSS.buildBalltree2(query, dimension, capacity, null, null, weight);
                datanode = indexDSS.buildBalltree2(dataMap, dimension, capacity, null, null, weight);
            }
        }

        Pair<Double, PriorityQueue<queueMain>> resultPair = AdvancedHausdorff.IncrementalDistance(query, dataMap, dimension, queryNode, datanode, 0, 1, 0, false, 0, false, queryindexmap, dataindexMap, null, true);
        long endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0 + "," + AdvancedHausdorff.disCompTime + "," + AdvancedHausdorff.EuclideanCounter);
        //	System.out.println(resultPair.getLeft());
        timeMultiple[0][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;

        // 2, top-k fast
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        HashMap<Integer, Double> result = Search.pruneByIndex(dataMapPorto, datasetRoot, queryNode, queryid,
                dimension, indexMap, datasetIndex, queryindexmap, datalakeIndex, datasetIdMapping, k,
                indexString, null, true, error, capacity, weight, saveDatasetIndex, query);
        endtime = System.nanoTime();
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;
    }

    /*
     * join and self-join
     */
    static void testJoinOutlier(int testOrder, int timeID, int queryID, int datasetID, Pair<Double, PriorityQueue<queueMain>> aPair) {
        numberofComparisons = 4;
        int countAlgorithm = 0;

        System.out.println("*******testing join and outlier");

        // 1. brute force, too slow
        long startTime1 = System.nanoTime();
        //	Join.scanning(dataMapPorto.get(queryID), dataMapPorto.get(datasetID), dimension);
        long endtime = System.nanoTime();
        System.out.println("join with top-1 dataset with scanning baseline costs: " + (endtime - startTime1) / 1000000000.0);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;

        // 2, using knn
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        Join.joinTableBaseline(dataMapPorto.get(queryID), dataMapPorto.get(datasetID),
                indexMap.get(queryID), indexMap.get(datasetID), dimension, indexDSS);
        endtime = System.nanoTime();
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;
        System.out.println("join with top-1 dataset with nn search baseline costs: " + (endtime - startTime1) / 1000000000.0);

        // 3, join with our cached queue
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        Join.IncrementalJoin(dataMapPorto.get(queryID), dataMapPorto.get(datasetID), dimension, indexMap.get(limit + 1),
                indexMap.get(datasetID), 0, 0, 0.01, false, 0, false, aPair.getLeft(), aPair.getRight(), null, null, "haus", null, true);
        endtime = System.nanoTime();
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;
        System.out.println("join with top-1 dataset on reused queues costs: " + (endtime - startTime1) / 1000000000.0);


        // 4, test self-join for outlier detection
        startTime1 = System.nanoTime();//exact search: 0, 1, 0
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        AdvancedHausdorff.setParameter(true, false);
        Pair<Double, PriorityQueue<queueMain>> resultPaira = AdvancedHausdorff.IncrementalDistance(dataMapPorto.get(1), dataMapPorto.get(1), dimension, indexMap.get(1), indexMap.get(1), 0, 1, 0, false, 0, false, null, null, null, true);
        endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0 + "," + AdvancedHausdorff.disCompTime + "," + AdvancedHausdorff.EuclideanCounter);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;
        //	System.out.println(resultPaira.getLeft());
        //	System.out.println();
    }


    // four methods: intersecting range, overlapped grids in scanning and pruning mode
//    测试range query
//    k: 返回多少个结果
//    testOrder：用来确定timeMultiple中的位置，暂时用不上
//    queryid：query的数据集的id，这个比较重要
    static void testRangeQuery(int k, int testOrder, int queryid, int timeID) throws FileNotFoundException, IOException {
        //scan every dataset to find the top-k
        numberofComparisons = 4;
        int countAlgorithm = 0;
        indexNode root = datasetRoot;//datalakeIndex.get(1);//use storee
        if (datalakeIndex != null)
            root = datalakeIndex.get(1);

        System.out.println("*******testing top-k range query");

        double[] querymax = new double[dimension];
        double[] querymin = new double[dimension];
        double[][] dataset = null;
        indexNode queryNode = getQueryNode(dataset, queryid);
        for (int i = 0; i < dimension; i++) {// to enlarge the range for evaluation
            querymax[i] = queryNode.getMBRmax()[i] + error;
            querymin[i] = queryNode.getMBRmin()[i] - error;
        }

        // 1, the query by intersecting area, scan every dataset
        Search.setScanning(true);
        HashMap<Integer, Double> result = new HashMap<>();
        long startTime1 = System.nanoTime();
        Search.rangeQueryRankingArea(root, result, querymax, querymin, Double.MAX_VALUE, k, null, dimension,
                datalakeIndex, dimNonSelected, dimensionAll);
        long endtime = System.nanoTime();
        System.out.println("range query cost time = " + (endtime - startTime1) / 1000000000.0);
//        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;

        /*// 2, the grid-based overlap
        int[] queryzcurve = new int[zcodemap.get(queryid).size()];
        startTime1 = System.nanoTime();
        for (int i = 0; i < zcodemap.get(queryid).size(); i++)
            queryzcurve[i] = zcodemap.get(queryid).get(i);
        EffectivenessStudy.topkAlgorithmZcurveHashMap(zcodemap, zcodemap.get(queryid), k);
        endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;

        // 3, use datalake index to rank by area
        startTime1 = System.nanoTime();
        Search.setScanning(false);
        result = new HashMap<>();
        Search.rangeQueryRankingArea(root, result, querymax,
                querymin, Double.MAX_VALUE, k, null, dimension, datalakeIndex, dimNonSelected, dimensionAll);
        endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;


        // 4, use datalake index to rank by number of overlapped grids
        startTime1 = System.nanoTime();
        result = new HashMap<>();
        Search.gridOverlap(root, result, queryzcurve, Double.MAX_VALUE, k, null, datalakeIndex);
        endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;*/
    }

    /*
     * we only test the varied resolution here
     */
    static void testGridOverlapAlone(int k, int testOrder, int queryid, int timeID) {
        numberofComparisons = 4;
        int countAlgorithm = 3;

        System.out.println("*******testing top-k grip overlap");

        int[] queryzcurve = new int[zcodemap.get(queryid).size()];
        long startTime1 = System.nanoTime();
        for (int i = 0; i < zcodemap.get(queryid).size(); i++)
            queryzcurve[i] = zcodemap.get(queryid).get(i);
        EffectivenessStudy.topkAlgorithmZcurveHashMap(zcodemap, zcodemap.get(queryid), k);
        long endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0);
//        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;
    }

    // range query by using dataset index, i.e., count how many points in the range, and plus return all the points
    static void testRangeQueryDeep(int k, int testOrder, int queryid, int timeID) throws FileNotFoundException, IOException {
        //scan every dataset to find the top-k
        numberofComparisons = 4;
        int countAlgorithm = 2;
        long startTime1 = System.nanoTime();

        System.out.println("*******testing top-k range query by covered points");

        double[] querymax = new double[dimension];
        double[] querymin = new double[dimension];
        for (int i = 0; i < dimension; i++) {// enlarge the range for evaluation
            querymax[i] = indexMap.get(queryid).getMBRmax()[i] + error;
            querymin[i] = indexMap.get(queryid).getMBRmin()[i] - error;
        }

        indexNode root = datasetRoot;//datalakeIndex.get(1);//use storee
        if (datalakeIndex != null)
            root = datalakeIndex.get(1);

        // 1, counting how many points in the range, access the dataset index
        ArrayList<double[]> points = new ArrayList<>();
        // find and store the founded points
        Search.setScanning(true);
        HashMap<Integer, Double> result = new HashMap<>();
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        Search.rangeQueryRankingNumberPoints(root, result, querymax,
                querymin, Double.MAX_VALUE, k, null, dimension, dataMapPorto, datalakeIndex,
                datasetIndex, datasetIdMapping, indexString, indexDSS, dimNonSelected, dimensionAll, points, indexMap, capacity, weight, saveDatasetIndex);
        long endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;

        // 2, stored the points in the range, access the dataset index
        points = null;
        startTime1 = System.nanoTime();
        Search.setScanning(false);
        result = new HashMap<>();
        startTime1 = System.nanoTime();// exact search: 0, 1, 0
        Search.rangeQueryRankingNumberPoints(root, result, querymax,
                querymin, Double.MAX_VALUE, k, null, dimension, dataMapPorto, datalakeIndex,
                datasetIndex, datasetIdMapping, indexString, indexDSS, dimNonSelected, dimensionAll, points, indexMap, capacity, weight, saveDatasetIndex);
        endtime = System.nanoTime();
        System.out.println((endtime - startTime1) / 1000000000.0);
        timeMultiple[timeID][testOrder * numberofComparisons + countAlgorithm++] += (endtime - startTime1) / 1000000000.0;
    }

    /*
     * running time record
     */
    static void recordRunningTime(String filename, int timeid) {
        String content = datalakeID + ",";
        for (double t : timeMultiple[timeid])
            content += t / NumberQueryDataset + ",";
        write(filename, content + "\n");
    }

    /*
     * running time record of
     */
    static void recordIndexModeEffect(String filename1, String filename2) {
        String content = datalakeID + ",";
        for (int i = 0; i < numberofParameterChoices; i++) {
            content += "0,";
            for (int j = 0; j < 3; j++) {
                double t = timeMultiple[j * 2][i * numberofComparisons + 3];
                content += t / NumberQueryDataset + ",";
            }
        }
        write(filename1, content + "\n");

        for (int i = 0; i < numberofParameterChoices; i++) {
            content += "0,";
            for (int j = 0; j < 3; j++) {
                double t = timeMultiple[j * 2 + 1][i * numberofComparisons + 3];
                content += t / NumberQueryDataset + ",";
            }
        }
        write(filename2, content + "\n");
    }

    /*
     * running time record
     */
    static void recordIndexPerformance(String filename, int timeid, double spaceMultiple[][]) {
        String content = datalakeID + ",";
        for (double t : spaceMultiple[timeid])
            content += t / NumberQueryDataset + ",";
        write(filename, content + "\n");
    }


    /*
     * update the datalake
     */
    static void addingDatasets() throws FileNotFoundException, IOException {
        /**
         * test adding new dataset into the datalake
         */
//		indexNodes = new ArrayList<indexNode>();
//		for(int dataid: datasetIndex.keySet()) {
//			indexNodes.add(datasetIndex.get(dataid).get(1));
//		}
//		indexDSS.updateDatasetNodeId(datasetRoot, 1);
        //	datalakeIndex = null;
        /* update index */
        //	indexNode rootBall = indexDSS.buildBalltree2(dataMapPorto.get(limit+1), dimension, capacity);
        //	rootBall.setroot(limit+1);
        //	datasetRoot = indexDSS.insertNewDataset(rootBall, datasetRoot, capacity, dimension);
        //	rootNum = indexDSS.getRootCount(datasetRoot);
        //	System.out.println("number of datasets: "+ rootNum);
        //	dataMapPorto = null;
        int k = 10;
        long startTime1 = System.nanoTime();
        // test top-k search, which also works for restored index.
        int queryid = 10;
        double[][] querydata = readSingleFile(datasetIdMapping.get(queryid));
        HashMap<Integer, Double> result = Search.pruneByIndex(dataMapPorto, datasetRoot, indexMap.get(queryid), queryid,
                dimension, indexMap, datasetIndex, datasetIndex.get(queryid), datalakeIndex, datasetIdMapping, k,
                indexString, null, true, error, capacity, weight, saveDatasetIndex, querydata);
        long endtime = System.nanoTime();
        System.out.println("top-1 search costs: " + (endtime - startTime1) / 1000000000.0);
        System.out.println(result + ", " + datasetIdMapping.get(result.entrySet().iterator().next().getKey()));
    }

    /*
     * pair-wise comparison, including Hausdorff, join, outlier, by varying the dataset scale, and recreate the index
     */
    static void EvaluatePairWise() throws IOException {
        System.out.println("evaluate pair-wise hausdorff");
        timeMultiple = new double[3][numberofComparisons * numberofParameterChoices];
        for (int i = 0; i < numberofParameterChoices; i++) { // test the scale
            for (int j = 0; j < NumberQueryDataset; j++) {
                int scale = ss[i];
                limit = 2 * scale;
                storeIndexMemory = false;
                readDatalake(limit);
                storeIndexMemory = true;
                indexNodes = new ArrayList<>();
                dataMapPorto = prepareLargeDatasets(scale); // generate datasets by combining more 10, 100, 1000, 10000 datasets
                long startTime1a = System.nanoTime();
                for (int a : dataMapPorto.keySet()) {//we can build index to accelerate
                    if (a > limit + NumberQueryDataset)
                        break;
                    createDatasetIndex(a, dataMapPorto.get(a)); // load index structure
                }
                long endtimea = System.nanoTime();
                write("./logs/spadas/efficiency/DatasetPairwiseIndexTime.txt", datalakeID + "," + limit + "," + dimension + "," + capacity + "," + (endtimea - startTime1a) / 1000000000.0 + "\n");
                write("./logs/spadas/efficiency/DatasetPairwiseIndexSize.txt", datalakeID + "," + limit + "," + dimension + "," + capacity + "," + (getDatalakeIndexSize(datasetRoot, false) + getAllDatasetIndexSize()) / (1024.0 * 1024) + "\n");

                // 1, the effect of dataset scale on pairwise
                Pair<Double, PriorityQueue<queueMain>> resultPair = testPairwiseHausdorff(1, 2, i, 0, false);// 1,

                // 2, the effect of dataset scale on joins
                testJoinOutlier(i, 1, 1, 2, resultPair); // 2,

                // 3, the effect of scale on range query
                error = spaceRange / ((int) Math.pow(2, rs[0]));
                testPairwiseRangeQuery(indexMap.get(1), dataMapPorto.get(1), i, 2);
            }
        }
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-n-haus.txt", 0);
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-n-join-outlier.txt", 1);
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-n-range-deep.txt", 1);

        System.out.println("evaluate pair-wise range");
        timeMultiple = new double[1][numberofComparisons * numberofParameterChoices];
        for (int i = 0; i < numberofParameterChoices; i++) { // test range
            for (int j = 0; j < NumberQueryDataset; j++) {
                resolution = rs[i];
                error = spaceRange / ((int) Math.pow(2, resolution));
                testPairwiseRangeQuery(indexMap.get(1), dataMapPorto.get(1), i, 0);
            }
        }
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-r-range-deep.txt", 0);
    }


    static void EvaluateHausCapacity() throws IOException {
        System.out.println("evaluate capacity's effect on hausdorff");
        timeMultiple = new double[1][numberofComparisons * numberofParameterChoices];
        for (int i = 0; i < numberofParameterChoices; i++) { // test the scale
            for (int j = 0; j < NumberQueryDataset; j++) {
                int scale = 100;
                if (datalakeID == 3)
                    scale = 1;
                limit = 2 * scale; // limit is changed here
                storeIndexMemory = false;
                readDatalake(limit);
                storeIndexMemory = true;
                indexNodes = new ArrayList<>();
                dataMapPorto = prepareLargeDatasets(scale); // generate datasets by combining more 10, 100, 1000, 10000 datasets
                capacity = fs[i];
                for (int a : dataMapPorto.keySet()) {//we can build index to accelerate
                    if (a > limit + NumberQueryDataset)
                        break;
                    createDatasetIndex(a, dataMapPorto.get(a)); // load index structure
                }
                testPairwiseHausdorff(1, 2, i, 0, true);// test fast only
            }
        }
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-f-haus.txt", 0);
    }


    /*
     *
     */
    private static void EvaluateDatalakeScale() throws FileNotFoundException, IOException {
        timeMultiple = new double[2][numberofComparisons * numberofParameterChoices];
        indexNodes = new ArrayList<>();
        System.out.println("evaluate datalake's scale");
        for (int Ni = 0; Ni < ns.length; Ni++) {// varying datalake scale, preparing datalake index, and evaluate top-k search
            int N = (int) (ns[Ni] * limit);
            int a = 1;
            for (indexNode node : indexNodesAll) {//
                if (a++ <= N) // we only index part of the dataset.
                    indexNodes.add(node);
            }

            long startTime1a = System.nanoTime();
            File tempFile = new File(indexString + "datalake" + N + "-" + capacity + "-" + dimension + ".txt");
            if (!tempFile.exists()) {
                storeIndexMemory = true;
            } else {
                storeIndexMemory = false;
            }
            createDatalake(N);// create or load the index of whole data lake, and store it into disk to avoid rebuilding it
            long endtimea = System.nanoTime();
            write("./logs/spadas/efficiency/DatalakeIndexTime.txt", datalakeID + "," + N + "," + dimension + "," + capacity + "," + (endtimea - startTime1a) / 1000000000.0 + "\n");
            write("./logs/spadas/efficiency/DatasetIndexSize.txt", datalakeID + "," + limit + "," + dimension + "," + capacity + "," + (getDatalakeIndexSize(datasetRoot, false) + getAllDatasetIndexSize()) / (1024.0 * 1024) + "\n");
            for (int j = 0; j < NumberQueryDataset; j++) {
                testTopkHausdorff(10, Ni, j + limit, 0, limit); // 6,
                testRangeQuery(10, Ni, j + limit, 1); // 7,
            }
        }
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-N-topk-haus.txt", 0);
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-N-range.txt", 1);
    }

    static void evaluatePairWiseHaus(int[][] ids) {
//        ids是一个m*2的二维数组，代表有m个测试用例，id[i][0]是ia，id[i][1]是ib
        double exactHDs = 0, oldAdvHDs = 0, newAdvHDs = 0, overlapHDs = 0;
        double exactTimes = 0, oldAdvTimes = 0, newAdvTimes = 0, overlapTimes = 0;
        for (int[] pair : ids) {
            int ia = pair[0];
            int ib = pair[1];

            System.out.println("ia = " + ia);
            System.out.println("ib = " + ib);

            indexNode aNode = indexNodes.get(ia);
            indexNode bNode = indexNodes.get(ib);
            double[][] aData = dataMapPorto.get(ia + 1);
            double[][] bData = dataMapPorto.get(ib + 1);

            long startTime;
            long endTime;
            double costTime;

//        List<Double> res = new ArrayList<>();

//        测试radius是否valid

//        exact hd
            startTime = System.nanoTime();
            double exactHD = Hausdorff.HausdorffExact(aData, bData, dimension, false);
            endTime = System.nanoTime();
            costTime = (endTime - startTime) / 1000000000.0;
            System.out.println("exact hd = " + exactHD);
            System.out.println("cost time = " + costTime);
            exactHDs += exactHD;
            exactTimes += costTime;

//        old advanced hd
            startTime = System.nanoTime();
            double oldAdvHD = AdvancedHausdorff.IncrementalDistance(aData, bData, dimension, aNode, bNode,
                    0, 1, 0, false, 0, true, null,
                    null, null, true).getLeft();
            endTime = System.nanoTime();
            costTime = (endTime - startTime) / 1000000000.0;
            System.out.println("old adv hd = " + oldAdvHD);
            System.out.println("cost time = " + costTime);
            oldAdvHDs += oldAdvHD;
            oldAdvTimes += costTime;

//        new advanced hd
            startTime = System.nanoTime();
            double newAdvHD = AdvancedHausdorff.advancedHausdorff(aData, bData, aNode, bNode, dimension).getLeft();
            endTime = System.nanoTime();
            costTime = (endTime - startTime) / 1000000000.0;
            System.out.println("new adv hd = " + newAdvHD);
            System.out.println("cost time = " + costTime);
            newAdvHDs += newAdvHD;
            newAdvTimes += costTime;

//            overlap hd
            startTime = System.nanoTime();
            OverlapHausdorff overlapHausdorff = new OverlapHausdorff();
            double overlapHD = overlapHausdorff.NOHD(aNode, bNode, aData, bData, dimension);
            endTime = System.nanoTime();
            costTime = (endTime - startTime) / 1000000000.0;
            System.out.println("overlap hd = " + overlapHD);
            System.out.println("cost time = " + costTime);
            overlapHDs += overlapHD;
            overlapTimes += costTime;
        }
        double newHDRatio = newAdvHDs / exactHDs;
        double newTimeRatio = newAdvTimes / exactTimes;
        double overlapHDRatio = overlapHDs / exactHDs;
        double overlapTimeRatio = overlapTimes / exactTimes;
        System.out.println("new hd ratio: " + newHDRatio);
        System.out.println("new time ration: " + newTimeRatio);
        System.out.println("overlap hd ratio: " + overlapHDRatio);
        System.out.println("overlap time ratio: " + overlapTimeRatio);
    }

    /*
     * evaluate the effect of f on the datalake index
     */
    static void EvaluateDatalakeCapacity() throws IOException {
        timeMultiple = new double[2][numberofComparisons * numberofParameterChoices];
        indexNodes = new ArrayList<>();
        int N = limit;
        int a = 1;
        for (indexNode node : indexNodesAll) {//
            if (a++ <= N) // we only index part of the dataset.
                indexNodes.add(node);
        }
        System.out.println("evaluate datalake index's capacity");
        for (int i = 0; i < numberofParameterChoices; i++) { // test the scale
            capacity = fs[i];
            long startTime1a = System.nanoTime();
            File tempFile = new File(indexString + "datalake" + N + "-" + capacity + "-" + dimension + ".txt");
            if (!tempFile.exists()) {
                storeIndexMemory = true;
            } else {
                storeIndexMemory = false;
            }
            createDatalake(N);// create or load the index of whole data lake, and store it into disk to avoid
            long endtimea = System.nanoTime();
            write("./logs/spadas/efficiency/DatalakeIndexTime.txt", datalakeID + "," + N + "," + dimension + "," + capacity + "," + (endtimea - startTime1a) / 1000000000.0 + "\n");
            write("./logs/spadas/efficiency/DatasetIndexSize.txt", datalakeID + "," + limit + "," + dimension + "," + capacity + "," + (getDatalakeIndexSize(datasetRoot, false) + getAllDatasetIndexSize()) / (1024.0 * 1024) + "\n");

            for (int j = 0; j < NumberQueryDataset; j++) {
                testTopkHausdorff(10, i, j + limit, 0, limit); // 6, test the top-k Hausdorff
                testRangeQuery(10, i, j + limit, 1); // 7, test the range query
            }
        }
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-f-topk.txt", 0);
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-f-range.txt", 1);
    }

    /*
     * evaluate the effect of dimension on pairwise Hausdorff
     */
    static void EvaluateHausDimension() throws IOException {
        System.out.println("evaluate dimenison's effect on hausdorff");
        if (datalakeID == 5) {
            timeMultiple = new double[1][numberofComparisons * numberofParameterChoices];
            for (int i = 0; i < numberofParameterChoices; i++) { // test the dimension
                for (int j = 0; j < NumberQueryDataset; j++) {
                    dimension = ds[i];
                    dimensionAll = false;
                    dimNonSelected = new boolean[11];
                    for (int p = dimension; p < 11; p++)
                        dimNonSelected[p] = true;// non selected
                    testPairwiseHausdorff(1, 2, i, 0, false);// test all
                }
            }
            recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-d-haus.txt", 0);
        }
    }

    /*
     * evaluate top-k range query which only use datalake
     * 4 comparisons, ranking by
     * intersecting area, and grid-based overlap
     */
//    用来测试range query的用时
//    干脆直接用testRangeQuery函数来测试好了！简单粗暴！
    static void EvaluateRangequery() throws FileNotFoundException, IOException {
        // range query by datalake index only
        System.out.println("evaluate range query");
        timeMultiple = new double[2][numberofComparisons * numberofParameterChoices];
        for (int i = 0; i < numberofParameterChoices; i++) { // test k
            for (int j = 0; j < NumberQueryDataset; j++) {
                // 1, the
                testRangeQuery(ks[i], i, j + limit, 0);
                //	testRangeQueryDeep(ks[i], i, j + limit, 1);
            }
        }
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-k-range.txt", 0);
        //	recordRunningTime("./logs/spadas/efficiency/"+datalakeID+"-k-range-deep.txt", 1);

        // range query by datalake index only
        timeMultiple = new double[2][numberofComparisons * numberofParameterChoices];
        for (int i = 0; i < numberofParameterChoices; i++) { // test range
            for (int j = 0; j < NumberQueryDataset; j++) {
                resolution = rs[i];
                error = spaceRange / ((int) Math.pow(2, resolution));
                testRangeQuery(10, i, j + limit, 0);
                //			testRangeQueryDeep(10, i, j + limit, 1);
            }
        }
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-r-range.txt", 0);
        //	recordRunningTime("./logs/spadas/efficiency/"+datalakeID+"-r-range-deep.txt", 1);

        if (datalakeID == 5) {
            timeMultiple = new double[2][numberofComparisons * numberofParameterChoices];
            for (int i = 0; i < numberofParameterChoices; i++) { // test only the range query, as z-curve does not work
                for (int j = 0; j < NumberQueryDataset; j++) {
                    dimension = ds[i];
                    dimensionAll = false;
                    dimNonSelected = new boolean[11];
                    for (int p = dimension; p < 11; p++)
                        dimNonSelected[p] = true;// non selected
                    testRangeQuery(10, i, j + limit, 0);
                    //				testRangeQueryDeep(10, i, j + limit, 1);
                }
            }
            recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-d-range.txt", 0);
            //		recordRunningTime("./logs/spadas/efficiency/"+datalakeID+"-d-range-deep.txt", 1);
            dimension = 11;
        }
    }

    /*
     * top-k hausdorff query by datasets
     */
    static void EvaluateTopkHaus() throws FileNotFoundException, IOException {
        System.out.println("evaluate top-k Haus");

        timeMultiple = new double[1][numberofComparisons * numberofParameterChoices];
        for (int i = 0; i < ks.length; i++) { // varying k
            for (int j = 0; j < NumberQueryDataset; j++)
                testTopkHausdorff(ks[i], i, j + limit, 0, limit);
        }
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-k-topk-haus.txt", 0);

        if (datalakeID == 5) {
            timeMultiple = new double[1][numberofComparisons * numberofParameterChoices];
            for (int i = 0; i < ds.length; i++) { // varying d
                for (int j = 0; j < NumberQueryDataset; j++) {
                    dimension = ds[i];
                    dimensionAll = false;
                    dimNonSelected = new boolean[11];
                    for (int p = dimension; p < 11; p++)
                        dimNonSelected[p] = true;// non selected
                    testTopkHausdorff(10, i, j + limit, 0, limit);
                }
            }
            recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-d-topk-haus.txt", 0);
        }

        timeMultiple = new double[1][numberofComparisons * numberofParameterChoices];
        for (int i = 0; i < numberofParameterChoices; i++) { // test range
            for (int j = 0; j < NumberQueryDataset; j++) {
                resolution = rs[i];
                error = spaceRange / ((int) Math.pow(2, resolution));
                testHausdorffApproxi(10, i, j + limit, 0);//cannot work for porto
            }
        }
        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-r-haus-appro.txt", 0);
    }

    private static void EvaluateGridResolution() throws FileNotFoundException, IOException {
        System.out.println("evaluate grid resolution on GBO query");
        dimension = 2;
        timeMultiple = new double[1][numberofComparisons * numberofParameterChoices];
        for (int Ri = 0; Ri < rs.length; Ri++) {
            resolution = rs[Ri];
            error = spaceRange / ((int) Math.pow(2, resolution));
            String zcurveFile = zcodeSer + resolution + "-" + limit + ".ser";
            File tempFile = new File(zcurveFile);
            if (!tempFile.exists()) {// create the z-curve code for each dataset
                zcurveExist = false;
                if (dataMapPorto == null)
                    zcodemap = EffectivenessStudy.storeZcurveFile(minx, miny, spaceRange, 2, resolution, zcurveFile, limit, datasetIdMapping);
                else
                    zcodemap = EffectivenessStudy.storeZcurve(minx, miny, spaceRange, 2, dataMapPorto, resolution, zcurveFile);
            } else {
                zcurveExist = true;
                zcodemap = EffectivenessStudy.deSerializationZcurve(zcurveFile);
            }
            if (datalakeIndex != null)
                datalakeIndex.get(1).buildsignarture(zcodemap, datalakeIndex);
            else {
                datasetRoot.buildsignarture(zcodemap, datalakeIndex);
            }
            for (int j = 0; j < NumberQueryDataset; j++) {
                testGridOverlapAlone(10, Ri, 100, 0);// 8, increase the grid-size to observe the performance
            }
        }
//        recordRunningTime("./logs/spadas/efficiency/" + datalakeID + "-r-topk-grid.txt", 0);
    }

    /*
     * evaluate the effect of index's memory mode on the performance
     */
    static void EvaluateIndex(int option) throws IOException {
        timeMultiple = new double[6][numberofComparisons * numberofParameterChoices];
        indexMultiple = new double[2][numberofComparisons * numberofParameterChoices];
        String fileNameString;
        if (option == 1) {
            fileNameString = "d";
        } else {
            fileNameString = "N";
        }
        for (int i = 0; i < numberofParameterChoices; i++) { // test the scale
            int N = limit; // or change to dimension test
            if (option == 1) {
                dimension = ds[i];
                dimensionAll = false;
                dimNonSelected = new boolean[11];
                for (int p = dimension; p < 11; p++)
                    dimNonSelected[p] = true;// non selected
            } else {
                N = (int) (ns[i] * limit) + 10;
            }
            System.out.println("datalake scale: " + N);
            long startTime1a = System.nanoTime();
            // ** run all in memory, and run top-k Hausdorff and range query
            storeAllDatasetMemory = true; // load most datasets and create index
            storeIndexMemory = true; // do not create index
            datalakeIndex = null;
            indexMap = null;
            dataMapPorto = null;
            readDatalake(N);// read data, create data set index
            if (dimension > 2)// read and set weight
                AdvancedHausdorff.setWeight(dimension, indexString + "weight.ser");
            indexNodesAll = new ArrayList<>(indexNodes);
            createDatalake(N); // create datalake index
            long endtimea = System.nanoTime();
            double size = (getDatalakeSize(dimension) + getDatalakeIndexSize(datasetRoot, false) + getAllDatasetIndexSize()) / (1024.0 * 1024.0);
            indexMultiple[0][i * numberofComparisons + 1] += size;
            System.out.println("1index size is: " + size);
            indexMultiple[1][i * numberofComparisons + 1] += (endtimea - startTime1a) / 1000000000.0;
            for (int j = 0; j < NumberQueryDataset; j++) {
                testTopkHausdorff(10, i, j + N, 0, N); // 6, test the top-k Hausdorff
                testRangeQuery(10, i, j + N, 1); // 7, test the range query
            }

            // store dataset in memory, load datalake index from disk
            startTime1a = System.nanoTime();
            readDatalake(N);
            File tempFile = new File(indexString + "datalake" + N + "-" + capacity + "-" + dimension + ".txt");
            if (!tempFile.exists()) {
                storeIndexMemory = true;
            } else {
                storeIndexMemory = false;
            }
            createDatalake(N);// create or load the index of whole data lake, and store it into disk to avoid
            endtimea = System.nanoTime();
            size = (getDatalakeSize(dimension) + getDatalakeIndexSizeDisk(dimension)) / (1024.0 * 1024);
            System.out.println("2index size is: " + size);
            indexMultiple[0][i * numberofComparisons + 2] += size;
            indexMultiple[1][i * numberofComparisons + 2] += (endtimea - startTime1a) / 1000000000.0;
            for (int j = 0; j < NumberQueryDataset; j++) {
                testTopkHausdorff(10, i, j + N, 2, N); // 6, test the top-k Hausdorff
                testRangeQuery(10, i, j + N, 3); // 7, test the range query
            }

            // just load the datalake index
            if (datalakeID <= 2)
                continue;
            dataMapPorto = null;
            datalakeIndex = null;
            storeAllDatasetMemory = false;
            buildOnlyRoots = true; // no dataset index
            startTime1a = System.nanoTime();
            tempFile = new File(indexString + "datalake" + N + "-" + capacity + "-" + dimension + ".txt");
            if (!tempFile.exists()) {
                storeIndexMemory = true;
            } else {
                storeIndexMemory = false;
            }
            createDatalake(N);// create or load the index of whole data lake, and store it into disk to avoid
            endtimea = System.nanoTime();
            System.out.println("size of index is " + datalakeIndex.size());
            size = (getDatalakeIndexSizeDisk(dimension)) / (1024.0 * 1024);
            System.out.println("3index size is: " + size);
            indexMultiple[0][i * numberofComparisons + 3] += size;
            indexMultiple[1][i * numberofComparisons + 3] += (endtimea - startTime1a) / 1000000000.0;
            for (int j = 0; j < NumberQueryDataset; j++) {
                testTopkHausdorff(10, i, j + N, 4, N); // 6, test the top-k Hausdorff
                testRangeQuery(10, i, j + N, 5); // 7, test the range query
            }
        }
        //write the log files, running time by extracting the columns
        recordIndexModeEffect("./logs/spadas/efficiency/" + datalakeID + "-M-" + fileNameString + "-topk.txt", "./logs/spadas/efficiency/" + datalakeID + "-M-" + fileNameString + "-range.txt");
        recordIndexPerformance("./logs/spadas/efficiency/" + datalakeID + "-" + fileNameString + "-MemorySize.txt", 0, indexMultiple);
        recordIndexPerformance("./logs/spadas/efficiency/" + datalakeID + "-" + fileNameString + "-ConstructionCost.txt", 1, indexMultiple);
    }

    /*
     * test all the functions related to Hausdorff, creating index,
     * pair-wise computation, top-k, join, outlier, restored index
     */
    public static void main(String[] args) throws IOException {

        //	generateIndexFeatures(); used for extract features
        if (args.length > 1)
            limit = Integer.valueOf(args[1]); //specifying the length
        int temp = limit;

        //	limit = 10000;
        System.out.println("testing datalake with scale " + limit);
//		指定要处理的数据集目录或文件
        aString = args[0]; // the dataset folder or file

        testFullDataset = true;

//		建立下层索引
        readDatalake(limit);
        testFullDataset = false;

//		建立上层索引
        createDatalake(limit);

//        评测range query
//        testRangeQuery(10, 0, 100, 0);

//        评测GBO top-k query
//        testGridOverlapAlone(10, 0, 100, 0);

//        AdvancedHausdorff.IncrementalDistance
//        测评Haus top-k query
//        int ia = 300, ib = 500;
//        indexNode aNode = indexNodes.get(ia);
//        indexNode bNode = indexNodes.get(ib);
//        double[][] aData = dataMapPorto.get(ia + 1);
//        double[][] bData = dataMapPorto.get(ib + 1);

        long startTime = System.nanoTime();
//        Pair<Double, PriorityQueue<queueMain>> advancedHDPair = AdvancedHausdorff.IncrementalDistance(aData, bData, dimension, aNode, bNode,
//                0, 1, 0, false, 0, true, null, null, null,
//                true);
//        double advancedHD = advancedHDPair.getLeft();
//        System.out.println("advanced HD = " + advancedHD);
        long endTime = System.nanoTime();
        double costTime = (endTime - startTime) / 1000000000.0;
//        System.out.println("cal advanced hd cost time : " + costTime);

//        Pair<Double, PriorityQueue<queueMain>> resultPair = AdvancedHausdorff.IncrementalDistance(aData, bData, dimension, aNode, bNode,
//                0, 1, 0, false, 0, true, null, null, null,
//                true);
//        startTime = System.nanoTime();// exact search: 0, 1, 0
//        Join.IncrementalJoin(aData, bData, dimension, indexMap.get(limit + 1),
//                indexMap.get(ia), 0, 0, 0.01, false, 0, false, resultPair.getLeft(), resultPair.getRight(), null, null, "haus", null, true);
//        endTime = System.nanoTime();
//
//        System.out.println("join with top-1 dataset on reused queues costs: " + (endTime - startTime) / 1000000000.0);
//
//        EvaluateGridResolution();
//
//        System.out.println("stop here!");


//		评测豪斯多夫距离计算
//		int ia = 449, ib = 519;
//		int ia = 217, ib = 525;
//		int ia = 217, ib = 295;
//		int ia = 217, ib = 235;
//		int ia = 0, ib = 4;

//		long startTime, endTime;
//		double costTime;
//
		OverlapHausdorff overlapHausdorff = new OverlapHausdorff();
		RandSampHausdorff randSampHausdorff = new RandSampHausdorff();
//		indexNode aNode = indexNodes.get(ia);
//		indexNode bNode = indexNodes.get(ib);
//		double[][] aData = dataMapPorto.get(ia + 1);
//		double[][] bData = dataMapPorto.get(ib + 1);
//
//		boolean isRadiusValid = testRadiusAll(aNode, aData);
//		if (!isRadiusValid) {
//			System.out.println("radius error!");
//		} else {
//			System.out.println("radius correct!");
//		}

        System.out.println();

//        测试有无outlier removal情况下，exact hd和new advcanced hd的准确度和用时的差别，用比值来表示
        Random rand = new Random();
        int[][] testIDs = new int[10][2];
        for (int i = 0; i < 10; i++) {
            int ia = rand.nextInt(indexNodes.size());
            int ib = rand.nextInt(indexNodes.size());
//            evaluatePairWiseHaus(ia, ib);
            testIDs[i][0] = ia;
            testIDs[i][1] = ib;
        }
//        evaluatePairWiseHaus(testIDs);

        testTopkHausdorff(10, 0, 217, 0, limit);
        testTopkHausdorff(10, 0, 217, 9, limit);

        System.out.println();

//        startTime = System.nanoTime();
////		calculate max radius
//        calculateMaxRadius();
//        endTime = System.nanoTime();
//        costTime = (endTime - startTime) / 1000000000.0;
//        System.out.println("calculateMaxRadius cost time: " + costTime);

//        test
        recordRadius("radius_before");
        startTime = System.nanoTime();
        double[][] data = calculateMaxRadiusByKneedleBefore();
        calculateMaxRadiusByKneedle();
        endTime = System.nanoTime();
        costTime = (endTime - startTime) / 1000000000.0;
        System.out.println("calculateMaxRadiusByKneedle cost time: " + costTime);

//        detect and delete outliers from index tree
//        startTime = System.nanoTime();
//        deleteOutlierFromAll();
//        endTime = System.nanoTime();
//        costTime = (endTime - startTime) / 1000000000.0;
//        System.out.println("deleteOutliers cost time: " + costTime);


//        update attributes
        startTime = System.nanoTime();
        updateIndexNodes();
        endTime = System.nanoTime();
        costTime = (endTime - startTime) / 1000000000.0;
        System.out.println("updateIndexNodes cost time: " + costTime);

        recordRadius("radius_after");
//        calculateMaxRadius();
        calculateMaxRadiusByKneedleAfter(data);
        calculateMaxRadiusByKneedle();
        endTime = System.nanoTime();
        deleteOutliersTime = (endTime - startTime) / 1000000000.0;
        System.out.println("delete outliers cost time: " + deleteOutliersTime);
        System.out.println("finish");

//        evaluatePairWiseHaus(testIDs);

        System.out.println();

        testTopkHausdorff(10, 0, 100, 0, limit);

        System.out.println();

//		overlapHausdorff.NOHD
//		startTime = System.nanoTime();
//		double overlapHD = overlapHausdorff.NOHD(aNode, bNode, aData, bData, dimension);
//		System.out.println("overlap hd = " + overlapHD);
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("cal overlap hd cost time : " + costTime);

//		Hausdorff.HausdorffExact
//		startTime = System.nanoTime();
//		double exactHD = Hausdorff.HausdorffExact(aData, bData, dimension, false);
//		System.out.println("exact HD = " + exactHD);
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("cal exact hd cost time : " + costTime);

//		Hausdorff.naiveHausdorff
//		startTime = System.nanoTime();
//		double naiveHD = Hausdorff.naiveHausdorff(aData, bData, dimension);
//		System.out.println("naive HD = " + naiveHD);
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("cal naive hd cost time : " + costTime);

//		AdvancedHausdorff.IncrementalDistance
//		startTime = System.nanoTime();
//		Pair<Double, PriorityQueue<queueMain>> advancedHDPair = AdvancedHausdorff.IncrementalDistance(aData, bData, dimension, aNode, bNode,
//				0, 1, 0, false, 0, true, null, null, null,
//				true);
//		double advancedHD = advancedHDPair.getLeft();
//		System.out.println("advanced HD = " + advancedHD);
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("cal advanced hd cost time : " + costTime);

//		Hausdorff.HausdorffWithIndexDFS
//		startTime = System.nanoTime();
//		double pivotDist = Util.EuclideanDis(aNode.getPivot(), bNode.getPivot(), dimension);
//		double ub = pivotDist + aNode.getRadius() + bNode.getRadius();
//		double knnHD = Hausdorff.HausdorffWithIndexDFS(aData, bData, dimension, aNode, bNode, 0, ub, indexDSS);
//		System.out.println("knn HD = " + knnHD);
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("cal knn hd cost time : " + costTime);

//		AdvancedHausdorff.advancedHausdorff
//		startTime = System.nanoTime();
//		Pair<Double, PriorityQueue<queueMain>> advancedHDPairNew = AdvancedHausdorff.advancedHausdorff(aData, bData, aNode, bNode,
//				dimension);
//		double advancedHDNew = advancedHDPairNew.getLeft();
//		System.out.println("advanced new HD = " + advancedHDNew);
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("cal advanced new hd cost time : " + costTime);

//		randSampHausdorff.run
//		startTime = System.nanoTime();
//		double randSampHD = randSampHausdorff.run(aNode, bNode, aData, bData, dimension);
//		System.out.println("random sample hd = " + randSampHD);
//		endTime = System.nanoTime();
//		costTime = (endTime - startTime) / 1000000000.0;
//		System.out.println("cal random sample hd cost time : " + costTime);


        // 1. pair-wise evaluation
        storeAllDatasetMemory = true; // load most datasets and create index
        storeIndexMemory = false; // do not create index
        zcurveExist = false;
        readDatalake(limit);//read data, create index if not
        if (dimension > 2)// read and set weight
            AdvancedHausdorff.setWeight(dimension, indexString + "weight.ser");
//        EvaluatePairWise(); // 1, preparing large datasets by increasing the scale, and evaluate pairwise comparison, slowest
//        EvaluateHausCapacity(); // 2, evaluate the capacity in pair-wise Hausdorff distance.


        limit = temp; // in case that limit changed
        //  2. top-k evaluation
        if (datalakeID <= 2)
            storeAllDatasetMemory = true;
        else {// datasets are large
            storeAllDatasetMemory = false; // no data in memory
            dataMapPorto = null;
            indexMap = null;
        }
        String zcurveFile = zcodeSer + resolution + "-" + limit + ".ser";
        File tempFile = new File(zcurveFile);
        if (!tempFile.exists()) {// create the z-curve code for each dataset
            zcurveExist = false;
        } else {
            zcurveExist = true;
            zcodemap = EffectivenessStudy.deSerializationZcurve(zcurveFile);
        }
        tempFile = new File(indexString + "datalake" + limit + "-" + capacity + "-" + dimension + ".txt");
        if (!tempFile.exists()) {
            storeIndexMemory = true;// has to create index as the datalake index is not there
            if (limit >= 100000 && datalakeID > 2)// for large datasets
                buildOnlyRoots = true;
        } else {
            storeIndexMemory = false;
            if (datalakeID <= 2)// has to create index for trajectory dataset
                storeIndexMemory = true;
        }
        long startTime1a = System.nanoTime(); // exact search: 0, 1, 0
        readDatalake(limit);// read the first
        indexNodesAll = new ArrayList<>(indexNodes);
        if (!zcurveExist)
            EffectivenessStudy.SerializedZcurve(zcurveFile, zcodemap);
        long endtimea = System.nanoTime();
//        write("./logs/spadas/efficiency/DatasetIndexTime.txt", datalakeID + "," + limit + "," + dimension + "," + capacity + "," + (endtimea - startTime1a) / 1000000000.0 + "\n");
        //	write("./logs/spadas/efficiency/DatasetIndexSize.txt", datalakeID+","+limit+","+dimension+","+capacity+","+(getAllDatasetIndexSize())/(1024.0*1024)+"\n");

        startTime1a = System.nanoTime();
        System.out.println(minx + "," + miny + "," + spaceRange + "," + resolution);
//        createDatalake(limit);// create or load the index of whole data lake, and store it into disk to avoid rebuilding it
        endtimea = System.nanoTime();
//        write("./logs/spadas/efficiency/DatalakeIndexTime.txt", datalakeID + "," + limit + "," + dimension + "," + capacity + "," + (endtimea - startTime1a) / 1000000000.0 + "\n");
//        write("./logs/spadas/efficiency/DatasetIndexSize.txt", datalakeID + "," + limit + "," + dimension + "," + capacity + "," + (getDatalakeIndexSize(datasetRoot, false) + getAllDatasetIndexSize()) / (1024.0 * 1024) + "\n");

//        EvaluateRangequery(); // 3, evaluate range query by increasing range
//        EvaluateTopkHaus(); // 4, evaluate top-k query by varying k, dimension, resolution
//        EvaluateHausDimension(); // 5, evaluate dimension

        if (datalakeID == 5) {// restore the parameter
            dimension = 11;
            dimensionAll = true;
        }

        // 3. varying datalake scale (number of datasets) to test
        limit = temp; // in case that limit changed
//        EvaluateDatalakeScale();

        // 4. capacity's effect to datalake index
        limit = temp; // in case that limit changed
//        EvaluateDatalakeCapacity();

        // 5. evaluate the resolution's effect on grid-based overlap
        limit = temp; // in case that limit changed
//        EvaluateGridResolution();


        // 6. evaluate the index mode's effect on construction and top-k Hausdorff and range query
        limit = temp;
//        EvaluateIndex(0); // changing N
//        if (datalakeID == 5)
//            EvaluateIndex(1); //*/// changing dimension
        System.out.println();
//		结束
    }

    /*
     * given a large dataset, simplified it into one dataset, and randomly sampled 99 datasets, integrate them and search by the large datasets.
     */
    void generateSimiplidedGroundTruth() {
        // use pick-means
    }
}
