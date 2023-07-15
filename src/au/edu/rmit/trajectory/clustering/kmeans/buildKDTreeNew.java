package au.edu.rmit.trajectory.clustering.kmeans;

import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.TimeUnit;

public class buildKDTreeNew {
    static double[][] dataset;
    static int[] index;
    static int minMergeLen = 5;
    static int maxStopLen = 15;
    static ForkJoinPool forkJoinPool;

    static void mergeSort(int[] index, int left, int right, int sk, int dim) {
        if (right - left > minMergeLen) {
            int mid = left + (right - left) / 2;
            mergeSort(index, left, mid, sk, dim);
            mergeSort(index, mid + 1, right, sk, dim);

            merge(index, left, mid, right, sk, dim);
        } else {
            insertSort(index, left, right, sk, dim);
//            System.out.println("nothing happened");
        }
    }

    static void insertSort(int[] index, int left, int right, int sk, int dim) {
        for (int i = left + 1; i < right + 1; i++) {
            int tmp = index[i];
            int j;
            for (j = i; j > left && compIndex(tmp, index[j - 1], sk, dim); j--) {
                index[j] = index[j - 1];
            }
            index[j] = tmp;
        }
    }

    static void merge(int[] index, int left, int mid, int right, int sk, int dim) {
        /*List<Integer> a = new ArrayList<>();
        List<Integer> b = new ArrayList<>();
        for (int i = left; i < mid + 1; i++) {
            a.add(index[i]);
        }
        a.add(Integer.MAX_VALUE);
        for (int j = mid + 1; j < right + 1; j++) {
            b.add(index[j]);
        }
        b.add(Integer.MAX_VALUE);

        for (int i = 0, j = 0, k = left; i < a.size() && j < b.size() && k < right + 1; k++) {
            index[k] = compIndex(a.get(i), b.get(j), sk, dim) ? a.get(i++) : b.get(j++);
        }*/

        int[] c = Arrays.copyOfRange(index, left, mid + 1);
        int[] d = Arrays.copyOfRange(index, mid + 1, right + 1);

        int i = 0, j = 0, k = left;
        for (; i < c.length && j < d.length; k++) {
            index[k] = compIndex(c[i], d[j], sk, dim) ? c[i++] : d[j++];
        }

        if (i == c.length) {
            while (j < d.length) {
                index[k++] = d[j++];
            }
        } else {
            while (i < c.length) {
                index[k++] = c[i++];
            }
        }
    }

    static boolean compIndex(int i1, int i2, int sk, int dim) {
        if (i1 == Integer.MAX_VALUE) {
            return false;
        } else if (i2 == Integer.MAX_VALUE) {
            return true;
        }

        double res = dataset[i1][sk] - dataset[i2][sk];

        for (int i = 1; res == 0 && i < dim; i++) {
            int r = (i + sk) % dim;
            res = dataset[i1][r] - dataset[i2][r];
        }

        return res < 0;
    }

    static KDTreeNode buildKDTree(int index[], int left, int right, int sk, int dim) {
        if (right - left > maxStopLen) {
            mergeSort(index, left, right, sk, dim);
            /*mergeTask task = new mergeTask(index, left, right, sk, dim);
            forkJoinPool.execute(task);
            try {
                forkJoinPool.awaitTermination(10, TimeUnit.MICROSECONDS);
            } catch (Exception e) {
                System.out.println(e.toString());
            }
*/
            int mid = left + (right - left) / 2;

            KDTreeNode root = new KDTreeNode(dataset, index, left, right, false);

            root.lChild = buildKDTree(index, left, mid, (sk + 1) % dim, dim);
            root.rChild = buildKDTree(index, mid + 1, right, (sk + 1) % dim, dim);

            return root;
        } else {
            return new KDTreeNode(dataset, index, left, right, true);
        }
    }

    static class mergeTask extends RecursiveAction {
        int maxTask = 100;
        int left;
        int right;
        int[] index;
        int sk;
        int dim;

        mergeTask(int[] index, int left, int right, int sk, int dim) {
            this.index = index;
            this.left = left;
            this.right = right;
            this.sk = sk;
            this.dim = dim;
        }

        protected void compute() {
            if (right - left < maxTask) {
                mergeSort(index, left, right, sk, dim);
            } else {
                int mid = left + (right - left) / 2;
                mergeTask leftTask = new mergeTask(index, left, mid, sk, dim);
                mergeTask rightTask = new mergeTask(index, mid + 1, right, sk, dim);

                leftTask.fork();
                rightTask.fork();

                leftTask.join();
                rightTask.join();

                merge(index, left, mid, right, sk, dim);
            }
        }
    }

    static indexNode createWonder(double[][] data, int dim) {
//        System.out.println("creat wonder!!!");
        dataset = data;
        index = new int[dataset.length];
        for (int i = 0; i < index.length; i++) {
            index[i] = i;
        }

        /*forkJoinPool = new ForkJoinPool();
        mergeTask task = new mergeTask(index, 0, index.length - 1, 0, dim);
        forkJoinPool.execute(task);
        try {
            forkJoinPool.awaitTermination(10, TimeUnit.MICROSECONDS);
        } catch (Exception e) {
            System.out.println(e.toString());
        }*/

        mergeSort(index, 0, index.length - 1, 0, dim);

        KDTreeNode root = buildKDTree(index, 0, index.length - 1, 0, dim);
        indexNode rootNode = new indexNode(dim);
        root.KDTreeToIndexTree(rootNode, dim);

        return rootNode;
    }

    public static void main(String[] args) throws InterruptedException {
//        double[][][] dataset = new double[2][20][3];
//        double[][] dataset = {{2,4,6},{4,1,6},{6,3,2},{4,6,5},{1,3,2},{7,3,3},{4,2,7},{8,2,3},{5,4,5},{4,1,7},
//                {2,5,3},{8,7,5},{7,8,4},{4,7,7},{2,5,1},{4,1,4},{3,7,4},{6,8,4},{6,5,3},{8,9,4},
//                {8,3,1},{7,4,6},{6,4,7},{7,4,3},{2,3,3},{5,3,7},{5,7,1},{3,5,6},{2,6,8},{3,3,3}};
        dataset = new double[100][3];
        Random rand = new Random(System.currentTimeMillis());
        for (int i = 0; i < dataset.length; i++) {
            for (int j = 0; j < dataset[0].length; j++) {
                dataset[i][j] = rand.nextInt(20);
            }
        }
        index = new int[100];
        /*for (int i = 0; i < index.length; i++) {
            index[i] = i;
        }

        ForkJoinPool forkJoinPool = new ForkJoinPool();
        mergeTask task = new mergeTask(index, 0, index.length - 1, 0, 3);
        forkJoinPool.execute(task);
        forkJoinPool.awaitTermination(2, TimeUnit.SECONDS);

//        buildKDTreeNew kdTreeNew = new buildKDTreeNew();
//        buildKDTreeNew.dataset = dataset;
//        mergeSort(index, 0, index.length - 1, 0, 3);
        KDTreeNode root = buildKDTree(index, 0, index.length - 1, 0, 3);
        indexNode rootNode = new indexNode(3);
        root.KDTreeToIndexTree(rootNode, 3);
        System.out.println("finish");*/
    }
}
