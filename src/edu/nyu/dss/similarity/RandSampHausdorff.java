package edu.nyu.dss.similarity;

import au.edu.rmit.trajectory.clustering.kmeans.indexNode;
import au.edu.rmit.trajectory.clustering.kpaths.Util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class RandSampHausdorff {
    int dimension;
    List<double[]> aData;
    List<double[]> bData;

    List<double[]> getDiffSet(indexNode bNode, double[][] a) {
        List<double[]> aNew = new ArrayList<>();
        double dist;
        double[] pivot = bNode.getPivot();
        double radius = bNode.getRadius();
        for (double[] p : a) {
            dist = Util.EuclideanDis(pivot, p, dimension);
            if (dist > radius) {
                aNew.add(p);
            }
        }
        return aNew;
    }

    void randomSample(indexNode bNode, double[][] a, double[][] b) {
        aData = getDiffSet(bNode, a);
        bData = new ArrayList<>();
        Collections.addAll(bData, b);
//        Collections.shuffle(aData);
//        Collections.shuffle(bData);
    }

    double run(indexNode aNode, indexNode bNode, double[][] a, double[][] b, int dimension) {
        this.dimension = dimension;
        randomSample(bNode, a, b);
        double maxDist = 0;
        for (double[] aPoint : aData) {
            double minDist = Double.MAX_VALUE;
            boolean isValid = true;
            for (double[] bPoint : bData) {
                double dist = Util.EuclideanDis(aPoint, bPoint, dimension);
                if (dist < maxDist) {
                    isValid = false;
                    break;
                }
                minDist = Math.min(minDist, dist);
            }
            if (isValid) {
                maxDist = Math.max(maxDist, minDist);
            }
        }
        return maxDist;
    }
}
