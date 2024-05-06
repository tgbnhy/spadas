package core.spadas.dss.similarity;

import tree.trajectory.clustering.kmeans.indexNode;

import java.util.*;

public class OverlapHausdorff {
    int dimension;
    double[][] A;
    double[][] B;

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

    double pointToPoint(double[] a, double[] b) {
        double dist = 0;
        for (int k = 0; k < dimension; k++) {
            dist += Math.pow((a[k] - b[k]), 2);
        }
        dist = Math.sqrt(dist);
        return dist;
    }

    double pointToPoints(double[] a, indexNode node, double maxLB) {
//        double[] a = A[i];
        double minDist = Double.MAX_VALUE;
        List<Integer> pointIdList = new ArrayList<>();
//        node.getPointIdListAll(pointIdList);
//        List<Integer> pointIdList = node.getpointIdList();
        for (double[] b : B) {
//            double[] b = B[j];
            double dist = pointToPoint(a, b);
            if (dist < maxLB) {
                minDist = dist;
                break;
            }
//            double dist = Util.EuclideanDis(a, b, dimension);
            minDist = Math.min(minDist, dist);
        }
        return minDist;
    }

    double lb(double[] p, indexNode node) {
        return Math.max(pointToPoint(p, node.getPivot()) - node.getRadius(), 0);
    }

    double ub(double[] p, indexNode node) {
        return pointToPoint(p, node.getPivot()) + node.getRadius();
    }

    double[] pointToCube(double[] p, indexNode node) {
        return new double[]{ub(p, node), lb(p, node)};
    }

//    aNode is cube/node, bNode is points/point set
    Map.Entry<double[], Boolean> cubeToPoints(indexNode aNode, indexNode bNode, double maxLB) {
        boolean isValid = true;
        double UB = Double.MAX_VALUE, LB = Double.MAX_VALUE;
        List<Integer> pointIdList = new ArrayList<>();
        bNode.getPointIdListAll(pointIdList);
        for (double[] b : B) {
//            double[] b = B[i];
            double ub = ub(b, aNode);
            double lb = lb(b, aNode);
            if (ub <= maxLB) {
                isValid = false;
                break;
            }
            UB = Math.min(UB, ub);
            LB = Math.min(LB, lb);
        }
        return new AbstractMap.SimpleEntry<>(new double[]{UB, LB}, isValid);
    }

    double NOHD(indexNode aNode, indexNode bNode, double[][] aData, double[][] bData, int dimension) {
        A = aData;
        B = bData;
        this.dimension = dimension;
//        Map.Entry<indexNode, Integer> entry = new AbstractMap.SimpleEntry<>();
//        descending queue
        Queue<Map.Entry<indexNode, Double>> dpq = new PriorityQueue<>((x1, x2) -> {
            return (int) (x2.getValue() - x1.getValue());
        });
        dpq.add(new AbstractMap.SimpleEntry<>(aNode, Double.MAX_VALUE));
        double maxLB = 0;
        while(!dpq.isEmpty()) {
            Map.Entry<indexNode, Double> entry = dpq.poll();
            indexNode node = entry.getKey();
            double minUB = entry.getValue();
            if (!node.isLeaf()) {
                if (minUB >= maxLB) {
                    Set<indexNode> nodeList = node.getNodelist();
                    for (indexNode cNode : nodeList) {
//                        cNode is cube(node), bNode is point set/points
                        Map.Entry<double[], Boolean> tmpEntry = cubeToPoints(cNode, bNode, maxLB);
                        double UB = tmpEntry.getKey()[0];
                        double LB = tmpEntry.getKey()[1];
                        boolean isValid = tmpEntry.getValue();
                        if (isValid) {
                            maxLB = Math.max(maxLB, LB);
                            dpq.add(new AbstractMap.SimpleEntry<>(cNode, UB));
                        }
                    }
                }
            } else {
                List<Integer> pointIdList = node.getpointIdList();
                for (int i : pointIdList) {
                    double[] a = A[i];
                    double dist = pointToPoints(a, bNode, maxLB);
                    if (dist < Double.MAX_VALUE) {
                        maxLB = Math.max(maxLB, dist);
                    }
                }
            }
        }
        return maxLB;
    }


}
