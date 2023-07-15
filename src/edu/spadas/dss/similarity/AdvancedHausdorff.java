package edu.spadas.dss.similarity;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.stream.Collectors;

import au.edu.rmit.trajectory.clustering.kpaths.Util;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;

import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import org.paukov.combinatorics3.Generator;

public class AdvancedHausdorff extends Hausdorff {

    static List<List<Integer>> permutations;

    static ArrayList<Integer> outlierID;// storing the outlier point id
    static ArrayList<Double> outlierDis;// storing the outlier point id
    static boolean selfOutlier = false; // only set as true when
    static int outlierThreshold = 0; // set as 0 for normal Hausdorff


    /* set the approximate version and outlier*/
    static void setParameter(boolean selfOutliers, boolean apporixmates) {
        selfOutlier = selfOutliers;
    }

    /*
     * based on radius, leaf (capacity), or depth, or covered points
     */
    static boolean stopSplitCondition(indexNode a, int option, int dimension, boolean nonselectedDimension[], boolean dimensionAll) {
        if (a == null)
            return true;

        if (option == 1) {
            if (a.getRadius(nonselectedDimension, dimension, dimensionAll, weight) > radiusThreshold)
                return false;
            else {
                return true;
            }
        } else if (option == 2) {
            if (a.getTotalCoveredPoints() > coveredPointsThreshold)
                return false;
            else {
                return true;
            }
        } else if (option == 3) {//this is for node b
            if (!a.isLeaf())
                return false;
            else {
                return true;
            }
        }
        return false;
    }

    /*
     * split the first node, or the second node based on
     */
    static boolean splitPriority(indexNode a, indexNode b, int option, boolean asplit, boolean bsplit, int dimension, boolean nonselectedDimension[], boolean dimensionAll) {
        if (asplit == true)
            return false;
        else if (bsplit == true)
            return true;

        if (option == 1) {//compute the similarity
            if (a.getRadius(nonselectedDimension, dimension, dimensionAll, weight) > b.getRadius(nonselectedDimension, dimension, dimensionAll, weight))
                return true;
            else {
                return false;
            }
        } else if (option == 2) {
            if (a.getTotalCoveredPoints() > b.getTotalCoveredPoints())
                return true;
            else {
                return false;
            }
        }
        return false;
    }

    /*
     * pami15, early breaking
     */
    public static double earlyBreaking(double[][] point1xys, double[][] point2xys, int dimension, boolean directed,
                                       ArrayList<Integer> aCoveredPoints, ArrayList<Integer> bCoveredPoints) {
        int x = aCoveredPoints.size();
        int y = bCoveredPoints.size();
        double[][] dist_matrix;
        dist_matrix = new double[x][y];
        double cmax = 0;
        for (int i = 0; i < x; i++) {
            double cmin = Double.MAX_VALUE;
            for (int j = 0; j < y; j++) {
                /*
                 * we add an early breaking on Euclidean distance, good for high-dimension dataset
                 */
                dist_matrix[i][j] = EuclideanDis(point1xys[aCoveredPoints.get(i) - 1], point2xys[bCoveredPoints.get(j) - 1], dimension);
                if (dist_matrix[i][j] < cmin)
                    cmin = dist_matrix[i][j];
                if (cmin < cmax)
                    break;
            }
            if (cmin > cmax && Double.MAX_VALUE > cmin)
                cmax = cmin;
        }

        if (directed == false)
            for (int j = 0; j < y; j++) {
                double cmin = Double.MAX_VALUE;
                for (int i = 0; i < x; i++) {
                    if (dist_matrix[i][j] == 0)
                        dist_matrix[i][j] = EuclideanDis(point1xys[aCoveredPoints.get(i) - 1], point2xys[bCoveredPoints.get(j) - 1], dimension);
                    if (dist_matrix[i][j] < cmin)
                        cmin = dist_matrix[i][j];
                    if (cmin < cmax)
                        break;
                }
                if (cmin > cmax && Double.MAX_VALUE > cmin)
                    cmax = cmin;
            }
        return cmax;
    }

    /*
     * using a brute force to compute the distance
     *
     * this is not efficient version, we need to have fast version without constructing the matrix
     */
    static double computeNodeDistance(indexNode a, indexNode b, int dimension, double[][] dataMatrixa, double[][] dataMatrixb,
                                      boolean direction, int apointid, int bpointid) {
        // get all the covered points of two nodes
        ArrayList<Integer> aCoveredPoints = new ArrayList<Integer>();
        if (apointid > 0)
            aCoveredPoints.add(apointid);
        else
            a.getAllCoveredPoints(aCoveredPoints);
        ArrayList<Integer> bCoveredPoints = new ArrayList<Integer>();
        if (bpointid > 0)
            bCoveredPoints.add(bpointid);
        else
            b.getAllCoveredPoints(bCoveredPoints);
        return earlyBreaking(dataMatrixa, dataMatrixb, dimension, direction, aCoveredPoints, bCoveredPoints);// take the early breaking idea
    }

    /*
     * our own method, with a general index, ball-tree, can work for any types of index tree.
     * An incremental Hausdorff distance calculation algorithm
     * 1) it can be used to detect own outliers, following the definition sigmod'00 for self outlier detection, we remove the self-pair
     * 2) support elbow method to return robust Hausdorff measure, leave it for hausdorff distance evaluation
     *
     * we need to input the first pair back to the queue
     */
    public static Pair<Double, PriorityQueue<queueMain>> IncrementalDistance(double[][] point1xys, double[][] point2xys,
                                                                             int dimension, indexNode X, indexNode Y, int splitOption, int fastMode,
                                                                             double error, boolean reverse, double directDis, boolean topkEarlyBreaking,
                                                                             Map<Integer, indexNode> nodelist, Map<Integer, indexNode> nodelist1,
                                                                             boolean nonselectedDimension[], boolean dimensionAll) {
        radiusThreshold = error;// the error for approximate search
        if (!reverse) {
            disCompTime = 0;
            EuclideanCounter = 0;
        }
        outlierID = new ArrayList<Integer>();
        outlierDis = new ArrayList<Double>();
        int counter = 0;
        if (splitOption == 0)
            fastMode = 0;//exact search, explore every point
        if (tightBound == 2)// MBR search
            permutations = Generator.permutation(0, 1).withRepetitions(dimension - 1).stream().collect(Collectors.<List<Integer>>toList());
        PriorityQueue<queueMain> aHeaps = new PriorityQueue<queueMain>();
        PriorityQueue<queueSecond> secondHeaps = new PriorityQueue<queueSecond>();
        queueSecond secondq = new queueSecond(Y, 0, 0);
        secondHeaps.add(secondq);
        double distance = 0;
        queueMain mainq = new queueMain(X, secondHeaps, Double.MAX_VALUE, 0);
        aHeaps.add(mainq);
        //the outlier can be put alone, here
        while (!aHeaps.isEmpty()) {
//    		mainq = aHeaps.poll();
            mainq = aHeaps.peek();
            secondHeaps = mainq.getQueue();
            secondq = secondHeaps.peek();// we need to think about here peek k continuous points, or we just set k=1
            indexNode anode = mainq.getIndexNode();
            double ub = mainq.getbound();
            if (reverse && ub <= directDis)// for the directed early breaking
                return new MutablePair<>(directDis, aHeaps);
            indexNode bnode = secondq.getNode();
            //	if(topkEarlyBreaking && anode!=null && bnode!=null && directDis<secondq.getbound())//for top-k search and pruning, estimate the bound
            //		continue;
            boolean asplit = stopSplitCondition(anode, splitOption);
            boolean bsplit = stopSplitCondition(bnode, splitOption);
            if (asplit == true && bsplit == true) {//if two nodes are null, i.e., they are points // we can change it to radius threshold
                if (topkEarlyBreaking) {// using radius to estimate the lower bound, if the lower bound still big, prune it,
                    double temp = ub - radiusThreshold;
                    if (anode == null && bnode == null)
                        return new MutablePair<>(ub, aHeaps);
                    if (bnode == null)
                        temp = ub;
                    if (directDis < temp) {// the estimated lower bound
                        return new MutablePair<>(ub, aHeaps);
                    } else {// otherwise, continue the exact search, change the split option to 0
                        splitOption = 0;
                        fastMode = 0;
                        if (anode != null)
                            asplit = false;
                        if (bnode != null)
                            bsplit = false;
                        topkEarlyBreaking = false;
                        if (asplit && bsplit == true) {
                            return new MutablePair<>(ub, aHeaps);
                        }
                    }
                } else {
                    // judge whether meet the stopping for outlier, and keep monitoring and infer the outlier
                    outlierID.add(mainq.getpointID());
                    outlierDis.add(ub);
                    // using the elbow method here to determine the precise value
                    if (++counter > outlierThreshold)
                        return new MutablePair<>(ub, aHeaps);
                }
            }
            if (asplit == false || bsplit == false) {
                boolean splitPr = splitPriority(anode, bnode, splitCompare, asplit, bsplit, dimension, nonselectedDimension, dimensionAll);//
                aHeaps.poll();
                if (splitPr) {
                    traverseX(anode, point1xys, point2xys, secondHeaps, dimension, aHeaps, fastMode, topkEarlyBreaking, directDis,
                            false, 0, nodelist, nonselectedDimension, dimensionAll);
                } else {
                    secondHeaps.poll();
                    ub = traverseY(bnode, point1xys, point2xys, secondHeaps, dimension, aHeaps, ub, anode, mainq.getpointID(),
                            fastMode, topkEarlyBreaking, directDis, false, 0, nodelist1, nonselectedDimension, dimensionAll);
                }
            }
        }
        Pair<Double, PriorityQueue<queueMain>> resultPair = new MutablePair<>(distance, aHeaps);
        return resultPair;
    }


    /*
     * reconstruct!!!!!
     * HausJoinAppr3In1!!!
     * viewed as a baseline
     */
//    噢原来这是我自己写的…

    public static Pair<Double, PriorityQueue<queueMain>> advancedHausdorff(double[][] qData, double[][] dData,
                                                                           indexNode qNode, indexNode dNode, int dimension) {
        PriorityQueue<queueSecond> pqSec = new PriorityQueue<>();
//		 maybe the third parameter pointId is useless
        pqSec.add(new queueSecond(dNode, 0, 0));
        PriorityQueue<queueMain> pqMain = new PriorityQueue<>();
//       maybe the third parameter pointId is useless
        pqMain.add(new queueMain(qNode, pqSec, Double.MAX_VALUE, 0));
        MutablePair<Double, PriorityQueue<queueMain>> res = new MutablePair<>();

        while (!pqMain.isEmpty()) {
            queueMain qMain = pqMain.peek();
            indexNode qNodeChild = qMain.getIndexNode();
            double ub = qMain.getbound();
            PriorityQueue<queueSecond> pqSecChild = qMain.getQueue();
            queueSecond qSecond = pqSecChild.peek();
            if (qSecond == null) {
                System.out.println();
            }
            indexNode dNodeChild = qSecond.getNode();
            if (qNodeChild != null && dNodeChild == null) {
                System.out.println();
            } else if (qNodeChild == null && dNodeChild != null){
                System.out.println();
            }
            if (qNodeChild == null && dNodeChild == null) {
//                double resDist = leafToLeaf(qData, dData, qNodeChild, dNodeChild, dimension);
                res.setLeft(ub);
                res.setRight(pqMain);
                break;
            } else {
                if (qNodeChild == null) {
                    pqMain.poll();
                    pqSecChild.poll();
                    traversalD(qNodeChild, dNodeChild, qData, dData, pqMain, pqSecChild, qMain, qSecond, ub, dimension);
                } else if (qNodeChild.getRadius() > dNodeChild.getRadius()) {
                    if (qNodeChild.isLeaf()) {
                        pqMain.poll();
                        pqSecChild.poll();
                        traversalD(qNodeChild, dNodeChild, qData, dData, pqMain, pqSecChild, qMain, qSecond, ub, dimension);
                    } else {
                        pqMain.poll();
                        traversalQ(qNodeChild, qData, dData, pqMain, pqSecChild, qMain, dimension);
                    }
                } else {
                    if (dNodeChild.isLeaf()) {
                        pqMain.poll();
                        traversalQ(qNodeChild, qData, dData, pqMain, pqSecChild, qMain, dimension);
                    } else {
                        pqMain.poll();
                        pqSecChild.poll();
                        traversalD(qNodeChild, dNodeChild, qData, dData, pqMain, pqSecChild, qMain, qSecond, ub, dimension);
                    }
                }
//                if ((!qNodeChild.isLeaf() && dNodeChild.isLeaf()) ||
//                        (!qNodeChild.isLeaf() && !dNodeChild.isLeaf() && qNodeChild.getRadius() > dNodeChild.getRadius())) {
//                    pqMain.poll();
//                    traversalQ(qNodeChild, dNodeChild, qData, dData, pqMain, qSecond, dimension);
//                } else {
//                    pqMain.poll();
//                    pqSecChild.poll();
//                    traversalD(qNodeChild, dNodeChild, qData, dData, pqMain, pqSecChild, ub, dimension);
//                }
//                if (!qNodeChild.isLeaf() && qNodeChild.getRadius() > dNodeChild.getRadius()) {
//                    pqMain.poll();
//                    traversalQ(qNodeChild, dNodeChild, qData, dData, pqMain, qSecond, dimension);
//                } else {
//                    pqMain.poll();
//                    pqSecChild.poll();
//                    traversalD(qNodeChild, dNodeChild, qData, dData, pqMain, pqSecChild, ub, dimension);
//                }
            }
        }
        return res;
    }

    public static double leafToLeaf(double[][] qData, double[][] dData, indexNode qNode, indexNode dNode, int dimension) {
        double cmax = 0;
        for (int i : qNode.getpointIdList()) {
            double cmin = Double.MAX_VALUE;
            for (int j : dNode.getpointIdList()) {
                double dist = Util.EuclideanDis(qData[i], dData[j], dimension);
                if (dist < cmax) {
                    cmin = 0;
                    break;
                }
                cmin = Math.min(cmin, dist);
            }
            cmax = Math.max(cmax, cmin);
        }
        return cmax;
    }

    //    q's radius > d's radius, so q cannot be a leaf. Just consider the situation that d is a leaf
    public static void traversalQ(indexNode qNode, double[][] qData, double[][] dData, PriorityQueue<queueMain> pqMain, PriorityQueue<queueSecond> pqSec, queueMain qMain, int dimension) {
//        List<queueSecond> listSec = new ArrayList<>();
//        while (!pqSec.isEmpty()) {
//            listSec.add(pqSec.poll());
//        }

        if (!qNode.isLeaf()) {
            for (indexNode qNodeChild : qNode.getNodelist()) {
                PriorityQueue<queueSecond> pqSecChild = new PriorityQueue<>();
                double minUB = Double.MAX_VALUE;
                for (queueSecond qSec : pqSec) {
                    indexNode dNode = qSec.getNode();
                    if (!dNode.isLeaf()) {
                        double[] bound = boundNodeToNode(qNodeChild, dNode, dimension);
                        double ub = bound[0];
                        double lb = bound[1];
                        pqSecChild.add(new queueSecond(dNode, lb, 0));
                        minUB = Math.min(minUB, ub);
                    } else {
                        double[] bound = boundNodeToLeaf(qNodeChild, dNode, dData, dimension);
                        double ub = bound[0];
                        double lb = bound[1];
                        pqSecChild.add(new queueSecond(dNode, lb, 0));
                        minUB = Math.min(minUB, ub);
                    }
                }
                pqMain.add(new queueMain(qNodeChild, pqSecChild, minUB, 0));
            }
        } else {
            for (queueSecond qSec : pqSec) {
                double minUB = Double.MAX_VALUE;
                indexNode dNode = qSec.getNode();
                PriorityQueue<queueSecond> pqSecChild = new PriorityQueue<>();
                if (dNode.isLeaf()) {
                    Pair<Double, int[]> pair = boundFromExact(qNode.getpointIdList(), dNode.getpointIdList(), qData, dData, dimension);
                    double dist = pair.getKey();
                    int[] pointID = pair.getValue();
                    pqSecChild.add(new queueSecond(null, dist, pointID[1]));
                    minUB = Math.min(minUB, dist);
                    pqMain.add(new queueMain(null, pqSecChild, minUB, pointID[0]));
                } else {
                    double[] bound = boundLeafToNode(qNode, dNode, qData, dimension);
                    double ub = bound[0];
                    double lb = bound[1];
                    pqSecChild.add(new queueSecond(dNode, lb, 0));
                    minUB = Math.min(minUB, ub);
                    pqMain.add(new queueMain(qNode, pqSecChild, minUB, 0));
                }
            }
        }

//        else {
//            double minUB = Double.MAX_VALUE;
//            for (queueSecond qSec : listSec) {
//                indexNode dNode = qSec.getNode();
//                Pair<Double, int[]> pair = new MutablePair<>();
//                PriorityQueue<queueSecond> pqSecChild = new PriorityQueue<>();
//                if (dNode == null) {
//                    pair = boundFromExact(qNode.getpointIdList(), new ArrayList<>(qSec.getPointId()), qData, dData, dimension);
//                    double dist = pair.getKey();
//                    int[] pointID = pair.getValue();
//                    pqSecChild.add(new queueSecond(null, dist, pointID[1]));
//                    minUB = Math.min(minUB, dist);
//                    pqMain.add(new queueMain(null, pqSecChild, minUB, pointID[0]));
//                } else if (!dNode.isLeaf()) {
//                    double[] bound = boundLeafToNode(qNode, dNode, qData, dimension);
//                    double ub = bound[0];
//                    double lb = bound[1];
//                    pqSecChild.add(new queueSecond(dNode, lb, 0));
//                    minUB = Math.min(minUB, ub);
//                    pqMain.add(new queueMain(qNode, pqSecChild, minUB, 0));
//                } else {
//                    pair = boundFromExact(qNode.getpointIdList(), dNode.getpointIdList(), qData, dData, dimension);
//                    double dist = pair.getKey();
//                    int[] pointID = pair.getValue();
//                    pqSecChild.add(new queueSecond(null, dist, pointID[1]));
//                    minUB = Math.min(minUB, dist);
//                    pqMain.add(new queueMain(null, pqSecChild, minUB, pointID[0]));
//                }
//            }
//        }


//        for (indexNode qNodeChild : qNode.getNodelist()) {
//            PriorityQueue<queueSecond> pqSecChild = new PriorityQueue<>();
//            double minUB = Double.MAX_VALUE;
//            if (dNode.isLeaf()) {
//                double lb = lbNodeToPoints(qNodeChild, dNode.getpointIdList(), dData, dimension);
//                double ub = ubNodeToPoints(qNodeChild, dNode.getpointIdList(), dData, dimension);
//                pqSecChild.add(new queueSecond(dNode, lb, 0));
//                minUB = Math.min(minUB, ub);
////                double lb = lbNodeToNode(qNodeChild, dNode, dimension);
////                double ub = ubNodeToNode(qNodeChild, dNode, dimension);
////                pqSecChild.add(new queueSecond(dNode, lb, 0));
////                minUB = Math.min(minUB, ub);
//            } else {
//                for (indexNode dNodeChild : dNode.getNodelist()) {
//                    double lb = lbNodeToNode(qNodeChild, dNodeChild, dimension);
//                    double ub = ubNodeToNode(qNodeChild, dNodeChild, dimension);
//                    pqSecChild.add(new queueSecond(dNodeChild, lb, 0));
//                    minUB = Math.min(minUB, ub);
//                }
//            }
//            pqMain.add(new queueMain(qNodeChild, pqSecChild, minUB, 0));
//        }
    }

    //    d's radius > q's radius, so d cannot be a leaf. Just consider the situation that q is a leaf
//    bullshit
    public static void traversalD(indexNode qNode, indexNode dNode, double[][] qData, double[][] dData, PriorityQueue<queueMain> pqMain, PriorityQueue<queueSecond> pqSec, queueMain qMain, queueSecond qSec, double UB, int dimension) {
        double minUB = UB;

        if (!dNode.isLeaf()) {
            if (qNode == null) {
                for (indexNode dNodeChild : dNode.getNodelist()) {
                    if (qMain.getpointID() < 0) {
                        System.out.println();
                    }
                    double[] bound = boundNullToNode(qMain.getpointID(), dNode, qData, dimension);
                    double ub = bound[0];
                    double lb = bound[1];
                    pqSec.add(new queueSecond(dNodeChild, lb, 0));
                    minUB = Math.min(minUB, ub);
                }
                pqMain.add(new queueMain(null, pqSec, minUB, qMain.getpointID()));
            } else {
                for (indexNode dNodeChild : dNode.getNodelist()) {
                    double[] bound = boundNodeToNode(qNode, dNodeChild, dimension);
                    double ub = bound[0];
                    double lb = bound[1];
                    pqSec.add(new queueSecond(dNodeChild, lb, 0));
                    minUB = Math.min(minUB, ub);
                }
                pqMain.add(new queueMain(qNode, pqSec, minUB, 0));
            }
        } else {
            Pair<Double, int[]> pair;
            if (qNode == null) {
                pair = boundFromExact(new ArrayList<>(qMain.getpointID()), new ArrayList<>(qSec.getPointId()), qData, dData, dimension);
                double dist = pair.getKey();
                int[] pointID = pair.getValue();
                pqSec.add(new queueSecond(null, dist, pointID[1]));
                minUB = Math.min(minUB, dist);
                pqMain.add(new queueMain(null, pqSec, minUB, pointID[0]));
            } else if (qNode.isLeaf()) {
                pair = boundFromExact(qNode.getpointIdList(), qSec.getNode().getpointIdList(), qData, dData, dimension);
                double dist = pair.getKey();
                int[] pointID = pair.getValue();
                pqSec.add(new queueSecond(null, dist, pointID[1]));
                minUB = Math.min(minUB, dist);
                pqMain.add(new queueMain(null, pqSec, minUB, pointID[0]));
            } else {
                double[] bound = boundNodeToLeaf(qNode, dNode, dData, dimension);
                double ub = bound[0];
                double lb = bound[1];
                pqSec.add(new queueSecond(dNode, lb, 0));
                minUB = Math.min(minUB, ub);
                pqMain.add(new queueMain(qNode, pqSec, minUB, 0));
            }
        }

//        if (dNode == null) {
//            Pair<Double, int[]> pair;
//            if (qNode == null) {
//                pair = boundFromExact(new ArrayList<>(qMain.getpointID()), new ArrayList<>(qSec.getPointId()), qData, dData, dimension);
//            } else {
//                pair = boundFromExact(qNode.getpointIdList(), new ArrayList<>(qSec.getPointId()), qData, dData, dimension);
//            }
//            double dist = pair.getKey();
//            int[] pointID = pair.getValue();
//            pqSec.add(new queueSecond(null, dist, pointID[1]));
//            minUB = Math.min(minUB, dist);
//            pqMain.add(new queueMain(null, pqSec, minUB, pointID[0]));
//        } else if (!dNode.isLeaf()) {
//            for (indexNode dNodeChild : dNode.getNodelist()) {
//                double[] bound = boundNodeToNode(qNode, dNodeChild, dimension);
//                double ub = bound[0];
//                double lb = bound[1];
//                pqSec.add(new queueSecond(dNodeChild, lb, 0));
//                minUB = Math.min(minUB, ub);
//            }
//            pqMain.add(new queueMain(qNode, pqSec, minUB, 0));
//        } else {
//            Pair<Double, int[]> pair;
//            if (qNode == null) {
//                pair = boundFromExact(new ArrayList<>(qMain.getpointID()), dNode.getpointIdList(), qData, dData, dimension);
//                double dist = pair.getKey();
//                int[] pointID = pair.getValue();
//                pqSec.add(new queueSecond(null, dist, pointID[1]));
//                minUB = Math.min(minUB, dist);
//                pqMain.add(new queueMain(null, pqSec, minUB, pointID[0]));
//            } else if (!qNode.isLeaf()) {
//                double[] bound = boundNodeToLeaf(qNode, dNode, dData, dimension);
//                double ub = bound[0];
//                double lb = bound[1];
//                pqSec.add(new queueSecond(dNode, lb, 0));
//                minUB = Math.min(minUB, ub);
//                pqMain.add(new queueMain(qNode, pqSec, minUB, 0));
//            } else {
//                if (!pqSec.isEmpty()) {
//                    pqMain.add(new queueMain(qNode, pqSec, minUB, 0));
//                }
//                pair = boundFromExact(qNode.getpointIdList(), dNode.getpointIdList(), qData, dData, dimension);
//                double dist = pair.getKey();
//                int[] pointID = pair.getValue();
//                PriorityQueue<queueSecond> pqSecNew = new PriorityQueue<>();
//                pqSecNew.add(new queueSecond(null, dist, pointID[1]));
//                minUB = Math.min(minUB, dist);
//                pqMain.add(new queueMain(null, pqSecNew, minUB, pointID[0]));
//            }
//        }

//        if (qNode.isLeaf()) {
//            for (indexNode dNodeChild : dNode.getNodelist()) {
//                double lb = lbPointsToNode(qNode.getpointIdList(), dNodeChild, qData, dimension);
//                double ub = ubPointsToNode(qNode.getpointIdList(), dNodeChild, qData, dimension);
//                pqSec.add(new queueSecond(dNodeChild, lb, 0));
//                minUB = Math.min(minUB, ub);
//            }
//            pqMain.add(new queueMain(qNode, pqSec, minUB, 0));
//        } else {
//            for (indexNode dNodeChild : dNode.getNodelist()) {
//                double lb = lbNodeToNode(qNode, dNodeChild, dimension);
//                double ub = ubNodeToNode(qNode, dNodeChild, dimension);
//                pqSec.add(new queueSecond(dNodeChild, lb, 0));
//                minUB = Math.min(minUB, ub);
//            }
//            pqMain.add(new queueMain(qNode, pqSec, minUB, 0));
//        }
//        for (indexNode dNodeChild : dNode.getNodelist()) {
//            double lb = lbNodeToNode(qNode, dNodeChild, dimension);
//            double ub = ubNodeToNode(qNode, dNodeChild, dimension);
//            pqSec.add(new queueSecond(dNodeChild, lb, 0));
//            minUB = Math.min(minUB, ub);
//        }
//        pqMain.add(new queueMain(qNode, pqSec, minUB, 0));
    }

    public static double lbNodeToPoints(indexNode a, List<Integer> pointIdList, double[][] dData, int dimension) {
        double lb = Double.MAX_VALUE;
        for (int i : pointIdList) {
            double temp = Math.max(Util.EuclideanDis(a.getPivot(), dData[i], dimension), 0);
            lb = Math.min(lb, temp);
        }
        return lb;
    }

    public static double ubNodeToPoints(indexNode a, List<Integer> pointIdList, double[][] dData, int dimension) {
        double ub = 0;
        for (int i : pointIdList) {
            double temp = Util.EuclideanDis(a.getPivot(), dData[i], dimension) + a.getRadius();
            ub = Math.max(ub, temp);
        }
        return ub;
    }

    public static double lbPointsToNode(List<Integer> pointIdList, indexNode b, double[][] qData, int dimension) {
        double lb = Double.MAX_VALUE;
        for (int i : pointIdList) {
            double temp = Math.max(Util.EuclideanDis(qData[i], b.getPivot(), dimension) - b.getRadius(), 0);
            lb = Math.min(lb, temp);
        }
        return lb;
    }

    public static double ubPointsToNode(List<Integer> pointIdList, indexNode b, double[][] qData, int dimension) {
        double ub = 0;
        for (int i : pointIdList) {
            double v1 = Util.EuclideanDis(qData[i], b.getPivot(), dimension);
            double temp = Math.sqrt(Math.pow(v1, 2) + Math.pow(b.getRadius(), 2));
            ub = Math.max(ub, temp);
        }
        return ub;
    }

    public static double lbNodeToNode(indexNode a, indexNode b, int dimension) {
        return Math.max(Util.EuclideanDis(a.getPivot(), b.getPivot(), dimension) - b.getRadius(), 0);
    }

    public static double ubNodeToNode(indexNode a, indexNode b, int dimension) {
        double v1 = Util.EuclideanDis(a.getPivot(), b.getPivot(), dimension);
        double v2 = Math.sqrt(Math.pow(v1, 2) + Math.pow(b.getRadius(), 2));
        return v2 + a.getRadius();
    }

    public static double[] boundNodeToNode(indexNode q, indexNode d, int dimension) {
        double[] res = new double[2];
        double v1 = Util.EuclideanDis(q.getPivot(), d.getPivot(), dimension);
        double v2 = Math.sqrt(Math.pow(v1, 2) + Math.pow(d.getRadius(), 2));
        res[0] = v2 + q.getRadius();
        res[1] = Math.max(Util.EuclideanDis(q.getPivot(), d.getPivot(), dimension) - d.getRadius(), 0);
        return res;
    }

    public static double[] boundNodeToLeaf(indexNode q, indexNode d, double[][] dData, int dimension) {
        double[] res = new double[2];
        double minDist = Double.MAX_VALUE;
        for (int i : d.getpointIdList()) {
            double[] p = dData[i];
            double tmp = Util.EuclideanDis(q.getPivot(), p, dimension);
            minDist = Math.min(minDist, tmp);
        }
        res[0] = minDist + q.getRadius();
        res[1] = Math.max(minDist - q.getRadius(), 0);
        return res;
    }

    public static double[] boundLeafToNode(indexNode q, indexNode d, double[][] qData, int dimension) {
        double[] res = new double[2];
        double maxDist = 0;
        for (int i : q.getpointIdList()) {
            double[] p = qData[i];
            double tmp = Util.EuclideanDis(p, d.getPivot(), dimension);
            maxDist = Math.max(maxDist, tmp);
        }
        res[0] = Math.sqrt(Math.pow(maxDist, 2) + Math.pow(d.getRadius(), 2));
        res[1] = Math.max(maxDist - d.getRadius(), 0);
        return res;
    }

    public static double[] boundNullToNode(Integer qID, indexNode dNode, double[][] qData, int dimension) {
        double[] res = new double[2];
        double v1 = Util.EuclideanDis(qData[qID], dNode.getPivot(), dimension);
        res[0] = Math.sqrt(Math.pow(v1, 2) + Math.pow(dNode.getRadius(), 2));
        res[1] = Math.max(v1 - dNode.getRadius(), 0);
        return res;
    }

    public static Pair<Double, int[]> boundFromExact(List<Integer> qPointIdList, List<Integer> dPointIdList, double[][] qData, double[][] dData, int dimension) {
//        dist, qPointID, dPointID
        double maxDist = 0;
        int maxQ = -1, minD = -1;
        for (int i : qPointIdList) {
            double minDist = Double.MAX_VALUE;
            for (int j : dPointIdList) {
                double dist = Util.EuclideanDis(qData[i], dData[j], dimension);
                if (dist < maxDist) {
                    minDist = 0;
                    break;
                }
                if (dist < minDist) {
                    minD = j;
                    minDist = dist;
                }
            }
            if (minDist > maxDist) {
                maxQ = i;
                maxDist = minDist;
            }
        }
        return new MutablePair<>(maxDist, new int[]{maxQ, minD});
    }

    /*
     * this is the directed version
     */
    static double IncrementalDistanceDirected(double[][] point1xys, double[][] point2xys, int dimension, indexNode X, indexNode Y,
                                              int splitOption, int fastMode, double error, Map<Integer, indexNode> nodelist, Map<Integer, indexNode> nodelist1, boolean nonselectedDimension[], boolean dimensionAll) {
        if (X.getRadius(nonselectedDimension, dimension, dimensionAll, weight) > Y.getRadius(nonselectedDimension, dimension, dimensionAll, weight)) {
            Pair<Double, PriorityQueue<queueMain>> resultPair = IncrementalDistance(point1xys, point2xys, dimension, X, Y,
                    splitOption, fastMode, error, false, 0, false, nodelist, nodelist1, nonselectedDimension, dimensionAll);
            Pair<Double, PriorityQueue<queueMain>> resultPairReverse = IncrementalDistance(point2xys, point1xys, dimension,
                    Y, X, splitOption, fastMode, error, true, resultPair.getLeft(), false, nodelist1, nodelist, nonselectedDimension, dimensionAll);
            return Math.max(resultPair.getLeft(), resultPairReverse.getLeft());
        } else {
            Pair<Double, PriorityQueue<queueMain>> resultPair = IncrementalDistance(point2xys, point1xys, dimension, Y, X,
                    splitOption, fastMode, error, false, 0, false, nodelist1, nodelist, nonselectedDimension, dimensionAll);
            Pair<Double, PriorityQueue<queueMain>> resultPairReverse = IncrementalDistance(point1xys, point2xys, dimension,
                    X, Y, splitOption, fastMode, error, true, resultPair.getLeft(), false, nodelist, nodelist1, nonselectedDimension, dimensionAll);
            return Math.max(resultPair.getLeft(), resultPairReverse.getLeft());
        }
    }

    // if it is outlier version, we exclude the self-pair
    static void traverseX(indexNode anode, double[][] point1xys, double[][] point2xys, PriorityQueue<queueSecond> secondHeaps,
                          int dimension, PriorityQueue<queueMain> aHeaps, int fastMode, boolean topkEarlyBreaking,
                          double directDis, boolean join, double joinThreshold, Map<Integer, indexNode> nodelist, boolean nonselectedDimension[], boolean dimensionAll) {
        ArrayList<double[]> pivtoList = new ArrayList<double[]>();
        ArrayList<Pair<double[], double[]>> mbrList = new ArrayList<Pair<double[], double[]>>();
        ArrayList<Double> radiusArrayList = new ArrayList<Double>();
        ArrayList<Integer> arrayList = new ArrayList<Integer>(anode.getpointIdList());
        ArrayList<indexNode> arrayListnode = new ArrayList<indexNode>(anode.getNodelist(nodelist));
        double newlb, newub;
        if (anode.isLeaf()) {
            for (int a : anode.getpointIdList()) {
//                原本是a - 1, 改成a
                pivtoList.add(point1xys[a]);
                radiusArrayList.add(0.0);
                if (tightBound == 2)
                    mbrList.add(new MutablePair<double[], double[]>(point1xys[a], point1xys[a]));
            }
        } else {
            for (indexNode a : anode.getNodelist(nodelist)) {
                pivtoList.add(a.getPivot());
                radiusArrayList.add(a.getRadius(nonselectedDimension, dimension, dimensionAll, weight));
                if (tightBound == 2)
                    mbrList.add(new MutablePair<double[], double[]>(a.getMBRmin(), a.getMBRmax()));
            }
        }
        int length = pivtoList.size();
        for (int i = 0; i < length; i++) {
            PriorityQueue<queueSecond> newsecondHeaps = new PriorityQueue<queueSecond>();
            double minub = Double.MAX_VALUE;
            for (queueSecond aQueueSecond : secondHeaps) {
                indexNode newnodeb = aQueueSecond.getNode();
                Map<Integer, Pair<double[], double[]>> segmentB = null;
                if (tightBound == 2) {//using the MBR to prune
                    if (newnodeb != null)
                        segmentB = generateSegment(newnodeb.getMBRmin(), newnodeb.getMBRmax(), dimension);
                    else
                        segmentB = generateSegment(point2xys[aQueueSecond.getPointId()], point2xys[aQueueSecond.getPointId()], dimension);
                }
                //if it has been computed and in the second round, we can retrieve the results and avoid the calculation
                double pivot_distance, bradius = 0;
                EuclideanCounter++;
                if (newnodeb != null) {
                    pivot_distance = newnodeb.getPivotdistance(pivtoList.get(i));
                    bradius = newnodeb.getRadius(nonselectedDimension, dimension, dimensionAll, weight);// can be changed to get radius
                } else {
//                    id - 1 => id
                    pivot_distance = Hausdorff.EuclideanDis(point2xys[aQueueSecond.getPointId()], pivtoList.get(i), dimension);
                }
                if (tightBound == 0) {
                    newlb = computeNewLowerBound(pivot_distance, radiusArrayList.get(i), bradius);
                    newub = computeNewUpperBound(pivot_distance, radiusArrayList.get(i), bradius);
                } else if (tightBound == 1) {
                    newlb = pivot_distance - bradius - radiusArrayList.get(i);
                    newub = pivot_distance + bradius + radiusArrayList.get(i);
                } else {
                    Pair<double[], double[]> mbrA = mbrList.get(i);//get the mbr box
                    Map<Integer, Pair<double[], double[]>> segmentA = generateSegment(mbrA.getLeft(), mbrA.getRight(), dimension);
                    newlb = computeMBRLowerBound(segmentA, segmentB, dimension);
                    newub = computeMBRUpperBound(segmentA, segmentB, dimension);
                    //	System.out.println(newlb+","+newub);
                }
                if (join && newlb > joinThreshold) {//if new lb is bigger than a given threshold, filter
                    continue;
                }

                if (fastMode == 1) {
                    if (newnodeb != null || !anode.isLeaf()) {
                        if (anode.isLeaf())
                            computeNodeDistance(null, newnodeb, dimension, point1xys, point2xys, true, arrayList.get(i), 0);
                        else
                            computeNodeDistance(arrayListnode.get(i), newnodeb, dimension, point1xys, point2xys, true, 0, aQueueSecond.getPointId());
                    }
                } else if (fastMode == 2) {
                    if (bradius < radiusThreshold && radiusArrayList.get(i) < radiusThreshold)
                        newub = pivot_distance; // approximate solution
                } else if (fastMode == 3) {
                    // one point and one leaf case, maybe we can compute in advance, as it is linear to compute rather than quadratic
                } else if (fastMode == 4) {
                    // learned distance based on embedding, very similar to RMI of the learned index, hybrid mode
                }

                if (!selfOutlier || arrayList.isEmpty() || arrayList.get(i) != aQueueSecond.getPointId()) {// self outlier detection, avoid self-pair, here
                    queueSecond secondqb = new queueSecond(newnodeb, newlb, aQueueSecond.getPointId());
                    newsecondHeaps.add(secondqb);
                    if (newub < minub)
                        minub = newub;
                }
            }
            if (anode.isLeaf()) {
                queueMain mainqa = new queueMain(null, newsecondHeaps, minub, arrayList.get(i));
                aHeaps.add(mainqa);
            } else {
                //	if(topkEarlyBreaking && newsecondHeaps.peek().getNode()!=null && directDis<newsecondHeaps.peek().getbound())//for top-k search and pruning, estimate the bound
                //		continue;
                queueMain mainqa = new queueMain(arrayListnode.get(i), newsecondHeaps, minub, 0);
                aHeaps.add(mainqa);
            }
        }
    }

    static double traverseY(indexNode bnode, double[][] point1xys, double[][] point2xys, PriorityQueue<queueSecond> secondHeaps,
                            int dimension, PriorityQueue<queueMain> aHeaps, double ub, indexNode anode, int apoint, int fastMode,
                            boolean topkEarlyBreaking, double directDis, boolean join, double joinThreshold, Map<Integer, indexNode> nodelist, boolean nonselectedDimension[], boolean dimensionAll) {
        ArrayList<double[]> pivtoList = new ArrayList<double[]>();
        ArrayList<Pair<double[], double[]>> mbrList = new ArrayList<Pair<double[], double[]>>();
        ArrayList<Integer> arrayList = new ArrayList<Integer>(bnode.getpointIdList());// this may cost much time
        ArrayList<indexNode> arrayListnode = new ArrayList<indexNode>(bnode.getNodelist(nodelist));
        ArrayList<Double> radiusArrayList = new ArrayList<Double>();
        double lb, newub;
        if (bnode.isLeaf()) {
            for (int a : bnode.getpointIdList()) {
//                a - 1 => a
                pivtoList.add(point2xys[a]);
                radiusArrayList.add(0.0);
                if (tightBound == 2)//mbr bound, complex
                    mbrList.add(new MutablePair<double[], double[]>(point2xys[a], point2xys[a]));
            }
        } else {
            for (indexNode a : bnode.getNodelist(nodelist)) {
                pivtoList.add(a.getPivot());
                radiusArrayList.add(a.getRadius(nonselectedDimension, dimension, dimensionAll, weight));
                if (tightBound == 2)
                    mbrList.add(new MutablePair<double[], double[]>(a.getMBRmin(), a.getMBRmax()));
            }
        }
        int length = pivtoList.size();
        Map<Integer, Pair<double[], double[]>> segmentA = null;
        if (tightBound == 2) {
            if (anode != null)
                segmentA = generateSegment(anode.getMBRmin(), anode.getMBRmax(), dimension);
            else
                //                apoint - 1 => apoint
                segmentA = generateSegment(point1xys[apoint], point1xys[apoint], dimension);
        }
        //	System.out.println(segmentA.size());
        for (int i = 0; i < length; i++) {
            double pivot_distance, bradius = 0;
            EuclideanCounter++;
            long startTime = System.nanoTime();
            if (anode != null) {
                pivot_distance = anode.getPivotdistance(pivtoList.get(i));
                bradius = anode.getRadius(nonselectedDimension, dimension, dimensionAll, weight);
            } else {
//                apoint - 1 => apoint
                pivot_distance = Hausdorff.EuclideanDis(point1xys[apoint], pivtoList.get(i), dimension);
            }
            long endtime = System.nanoTime();
            disCompTime += (endtime - startTime) / 1000000000.0;
            if (tightBound == 0) {
                lb = computeNewLowerBound(pivot_distance, bradius, radiusArrayList.get(i));
                newub = computeNewUpperBound(pivot_distance, bradius, radiusArrayList.get(i));
            } else if (tightBound == 1) {
                lb = pivot_distance - bradius - radiusArrayList.get(i);
                newub = pivot_distance + bradius + radiusArrayList.get(i);
            } else {
                Pair<double[], double[]> mbrB = mbrList.get(i);//get the mbr box
                Map<Integer, Pair<double[], double[]>> segmentB = generateSegment(mbrB.getLeft(), mbrB.getRight(), dimension);
                lb = computeMBRLowerBound(segmentA, segmentB, dimension);
                newub = computeMBRUpperBound(segmentA, segmentB, dimension);
            }
            if (join && lb > joinThreshold)// used for join only
                continue;
            boolean self = false;
            if (bnode.isLeaf()) {
                queueSecond mainqa = new queueSecond(null, lb, arrayList.get(i));
                if (!selfOutlier || arrayList.get(i) != apoint)// self outlier detection, avoid self-pair points
                    secondHeaps.add(mainqa);
                else {
                    self = true;
                }
            } else {
                queueSecond mainqa = new queueSecond(arrayListnode.get(i), lb, 0);
                secondHeaps.add(mainqa);
            }
            if (fastMode == 1) {//only when early access is used
                if (anode != null || !bnode.isLeaf()) {//no need to recompute if two are points
                    if (bnode.isLeaf())
                        computeNodeDistance(anode, null, dimension, point1xys, point2xys, true, 0, arrayList.get(i));
                    else
                        computeNodeDistance(anode, arrayListnode.get(i), dimension, point1xys, point2xys, true, apoint, 0);
                }
            } else if (fastMode == 2) {//the approximate version
                if (bradius < radiusThreshold && radiusArrayList.get(i) < radiusThreshold)
                    newub = pivot_distance; // for generating ground truth in a fast way, and estimating the lower bound based on the predicted distance
            } else if (fastMode == 3) {
                // one point and one leaf case, maybe we can compute in advance, as it is linear to compute rather than quadratic
            } else if (fastMode == 4) {
                // learned distance based on embedding, very similar to RMI (learned index) of the learned index, hybrid mode
            }
            if (!self && newub < ub)
                ub = newub;
        }
        //	if(topkEarlyBreaking && anode!= null && secondHeaps.peek().getNode()!=null && directDis<secondHeaps.peek().getbound())//for top-k search and pruning, estimate the bound
        //		return ub;
        queueMain mainqa = new queueMain(anode, secondHeaps, ub, apoint);
        aHeaps.add(mainqa);
        return ub;
    }

    /*
     * we propose tight upper bounds based half-ball theory, which is tighter than using whole-ball
     */
    static double computeNewUpperBound(double pivotDis, double radiusA, double radiusB) {
        return Math.sqrt((pivotDis) * (pivotDis) + radiusB * radiusB) + radiusA;
    }

    /*
     * we use tighter lower bounds based on balls
     */
    static double computeNewLowerBound(double pivotDis, double radiusA, double radiusB) {
        return pivotDis - radiusB;
    }

    /*
     * we compute the MBR upper bound, for general cases, i.e., any dimensional data
     */
    static double computeMBRUpperBound(Map<Integer, Pair<double[], double[]>> segmentA, Map<Integer, Pair<double[], double[]>> segmentB, int dimension) {
        double max = 0;
        for (int segid : segmentA.keySet()) {
            Pair<double[], double[]> aaPair = segmentA.get(segid);
            double min = Double.MAX_VALUE;
            for (int segidB : segmentB.keySet()) {
                Pair<double[], double[]> bbPair = segmentB.get(segidB);
                double mindist = segmentMaxDis(aaPair, bbPair, dimension);
                if (mindist < min)
                    min = mindist;
            }
            if (min > max)
                max = min;
        }
        return max;
    }

    /*
     * we compute the MBR lower bound, for general cases, i.e., any dimensional data
     */
    static double computeMBRLowerBound(Map<Integer, Pair<double[], double[]>> segmentA, Map<Integer, Pair<double[], double[]>> segmentB, int dimension) {
        double max = 0;
        for (int segid : segmentA.keySet()) {
            Pair<double[], double[]> aaPair = segmentA.get(segid);
            double min = Double.MAX_VALUE;
            for (int segidB : segmentB.keySet()) {
                Pair<double[], double[]> bbPair = segmentB.get(segidB);
                double mindist = segmentMinDis(aaPair, bbPair, dimension);
                if (mindist < min)
                    min = mindist;
            }
            if (min > max)
                max = min;
        }
        return max;
    }

    /*
     * generate all the edges of an MBR in d-dimensions
     */
    static Map<Integer, Pair<double[], double[]>> generateSegment(double[] mbrMinA, double[] mbrMaxA, int dimension) {
        Map<Integer, Pair<double[], double[]>> segmentMapA = new HashMap<Integer, Pair<double[], double[]>>();
        int segcount = 0;
        for (int i = 0; i < dimension; i++) {
            double min = mbrMinA[i];
            double max = mbrMaxA[i];
            for (List<Integer> point : permutations) {
                double[] pointA = new double[dimension];
                double[] pointB = new double[dimension];
                for (int j = 0; j < dimension; j++) {
                    if (j > i) {
                        if (point.get(j - 1) == 0) {
                            pointA[j] = mbrMaxA[j];
                            pointB[j] = mbrMaxA[j];
                        } else {
                            pointA[j] = mbrMinA[j];
                            pointB[j] = mbrMinA[j];
                        }
                    } else if (j == i) {
                        pointA[j] = min;
                        pointB[j] = max;
                    } else {
                        if (point.get(j) == 0) {
                            pointA[j] = mbrMaxA[j];
                            pointB[j] = mbrMaxA[j];
                        } else {
                            pointA[j] = mbrMinA[j];
                            pointB[j] = mbrMinA[j];
                        }
                    }
                }
                segmentMapA.put(segcount++, new MutablePair<>(pointA, pointB));
            }
        }
        return segmentMapA;
    }

    /*
     * compute the minimum distance between two segments by finding nearest of two points, from each other
     */
    static double segmentMinDis(Pair<double[], double[]> aaPair, Pair<double[], double[]> bbPair, int dimension) {
        double a = pointToSegmentDis(aaPair.getLeft(), bbPair, dimension);
        double b = pointToSegmentDis(aaPair.getRight(), bbPair, dimension);
        double c = pointToSegmentDis(bbPair.getLeft(), aaPair, dimension);
        double d = pointToSegmentDis(bbPair.getRight(), aaPair, dimension);
        return Math.min(a, Math.min(b, Math.min(c, d)));
    }

    /*
     * compute the minimum distance between a point and a segment
     */
    static double pointToSegmentDis(double[] start, Pair<double[], double[]> bbPair, int dimension) {
        double[] point1 = bbPair.getLeft();
        double[] point2 = bbPair.getRight();
        int di = 0;
        for (int d = 0; d < dimension; d++) {
            if (point1[d] != point2[d]) {
                di = d;
                break;
            }
        }
        if (start[di] >= Math.min(point1[di], point2[di]) && start[di] <= Math.max(point1[di], point2[di])) {
            point1[di] = start[di];
            return Hausdorff.EuclideanDis(start, point1, dimension);
        } else {
            return Math.min(Hausdorff.EuclideanDis(start, point1, dimension), Hausdorff.EuclideanDis(start, point2, dimension));
        }
    }

    /*
     * compute the maximum distance between two segments, can be faster
     */
    static double segmentMaxDis(Pair<double[], double[]> aaPair, Pair<double[], double[]> bbPair, int dimension) {
        double a = Hausdorff.EuclideanDis(aaPair.getLeft(), bbPair.getLeft(), dimension);
        double b = Hausdorff.EuclideanDis(aaPair.getRight(), bbPair.getRight(), dimension);
        double c = Hausdorff.EuclideanDis(aaPair.getLeft(), bbPair.getRight(), dimension);
        double d = Hausdorff.EuclideanDis(aaPair.getRight(), bbPair.getLeft(), dimension);
        return Math.max(a, Math.max(b, Math.max(c, d)));
    }

}
