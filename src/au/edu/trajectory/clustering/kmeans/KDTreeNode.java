package au.edu.rmit.trajectory.clustering.kmeans;

import au.edu.rmit.trajectory.clustering.kpaths.Util;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class KDTreeNode {
    Set<Integer> coveredPoints;
    double[][] dataset;
    int[] index;
    double[] pivot;
    double radius;
    double[] sum;
    boolean isLeaf;
    KDTreeNode lChild;
    KDTreeNode rChild;

    KDTreeNode(double[][] dataset, int[] index, int left, int right, boolean isLeaf) {
        this.dataset = dataset;
        this.index = index;
        this.isLeaf = isLeaf;
        pivot = new double[dataset[0].length];
        sum = new double[dataset[0].length];

        coveredPoints = new HashSet<>();
        for (int i = left; i < right + 1; i++) {
            coveredPoints.add(index[i]);
            for (int j = 0; j < sum.length; j++) {
                sum[j] += dataset[index[i]][j];
            }
        }
        for (int i = 0; i < pivot.length; i++) {
            pivot[i] = sum[i] / coveredPoints.size();
        }

        radius = 0;
        for (int p : coveredPoints) {
            radius = Math.max(radius, Util.EuclideanDis(pivot, dataset[p], pivot.length));
        }
    }

    void KDTreeToIndexTree(indexNode root, int dim) {
        root.setPivot(pivot);
        root.setRadius(radius);
        root.setSum(sum);
        root.setTotalCoveredPoints(coveredPoints.size());
        if (isLeaf) {
            root.addPoint(coveredPoints);
        } else {
            indexNode leftNode = new indexNode(dim);
            root.addNodes(leftNode);
            lChild.KDTreeToIndexTree(leftNode, dim);

            indexNode rightNode = new indexNode(dim);
            root.addNodes(rightNode);
            rChild.KDTreeToIndexTree(rightNode, dim);
        }
    }
}
