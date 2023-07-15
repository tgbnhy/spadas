package core.spadas.dss.similarity;

import tree.trajectory.clustering.kmeans.indexNode;
import tree.trajectory.clustering.kpaths.Util;

import java.util.*;

public class RandSampHausdorffOld {
    int dimension;
    double[][] A;
    double[][] B;

    List<indexNode> aLeafNodes;
    List<indexNode> bLeafNodes;

    double pointToPoint(int i, int j) {
        double[] a = A[i];
        double[] b = B[j];
        double dist = 0;
        for (int k = 0; k < dimension; k++) {
            dist += Math.pow((a[k]-b[k]), 2);
        }
        dist = Math.sqrt(dist);
        return dist;
    }

    double NodeToNode(List<Integer> aList, List<Integer> bList, double maxDist) {
        double hd = 0;
        for (int i : aList) {
            double minDist = Double.MAX_VALUE;
            for (int j : bList) {
                double dist = pointToPoint(i, j);
                if (dist < maxDist) {
                    break;
                }
                minDist = Math.min(minDist, dist);
            }
            if (minDist < Double.MAX_VALUE) {
                hd = Math.max(hd, minDist);
            }
        }
        if (hd == Double.MAX_VALUE) {
            System.out.println("something wrong!");
        }
        return hd;


    }

    void getLeafNodes(indexNode node, List<indexNode> leafNodes) {
        if (node.isLeaf()) {
            leafNodes.add(node);
        } else {
            for (indexNode child : node.getNodelist()) {
                getLeafNodes(child, leafNodes);
            }
        }
    }

    List<indexNode> randomSample(indexNode node) {
        List<indexNode> leafNodes = new ArrayList<>();
        getLeafNodes(node, leafNodes);
        Collections.shuffle(leafNodes);
        return leafNodes;
    }

    double getLowerBound(indexNode aNode, indexNode bNode) {
        return Math.max(Util.EuclideanDis(aNode.getPivot(), bNode.getPivot(), dimension) - bNode.getRadius(), 0);
    }

    double getUpperBound(indexNode aNode, indexNode bNode) {
        double v1 = Util.EuclideanDis(aNode.getPivot(), bNode.getPivot(), dimension);
        double v2 = Math.sqrt(Math.pow(v1, 2) + Math.pow(bNode.getRadius(), 2));
        return v2 + aNode.getRadius();
    }

    double run(indexNode aNode, indexNode bNode, double[][] aData, double[][] bData, int dimension) {
        A = aData;
        B = bData;
        this.dimension = dimension;
        aLeafNodes = randomSample(aNode);
        bLeafNodes = randomSample(bNode);
        double maxDist = 0;
        for (indexNode aChild : aLeafNodes) {
            double minDist = Double.MAX_VALUE;
            for (indexNode bChild : bLeafNodes) {
                double ub = getUpperBound(aChild, bChild);
                if (ub < maxDist) {
                    break;
                }
                double lb = getLowerBound(aChild, bChild);
                if (lb < minDist) {
//                    change to calculate the min min rather than the max min
                    double hd = NodeToNode(aChild.getpointIdList(), bChild.getpointIdList(), maxDist);
                    minDist = Math.min(minDist, hd);
                }
            }
            if (minDist < Double.MAX_VALUE) {
                maxDist = Math.max(maxDist, minDist);
            }
        }
        return maxDist;
    }
}
