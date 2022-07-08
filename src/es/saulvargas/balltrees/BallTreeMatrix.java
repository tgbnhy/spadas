/* 
 * Copyright (C) 2015 Saúl Vargas http://saulvargas.es
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package es.saulvargas.balltrees;

import it.unimi.dsi.fastutil.ints.IntArrayList;

import static java.lang.Math.max;
import static java.lang.Math.sqrt;

<<<<<<< HEAD
import java.util.Arrays;
=======
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
import java.util.HashSet;
import java.util.Random;
import java.util.Set;


import au.edu.rmit.trajectory.clustering.kmeans.indexNode;

/**
 * Ball tree.
 *
 * @author Saúl Vargas (Saul.Vargas@glasgow.ac.uk)
 */
public class BallTreeMatrix extends BinaryTree {

    private static final Random random = new Random();
    static double weight[];// indicate the weight on each dimension, for normalization
    
    public BallTreeMatrix(NodeBall root) {
        super(root);
    }
    
    public static void setWeight(int dimension, double []weightinput) {
    	weight = new double[dimension];
		for(int i=0; i< dimension; i++)
			if(weightinput==null)
				weight[i] = 1.0;
			else
				weight[i] = weightinput[i];
    }

    @Override
    public Ball getRoot() {
        return (Ball) super.getRoot();
    }

    public static indexNode create(double[][] itemMatrix, int leafThreshold, int maxDepth) {
        int[] rows = new int[itemMatrix.length];
        for (int row = 0; row < itemMatrix.length; row++) {
            rows[row] = row;
        }
        Ball root = new Ball(rows, itemMatrix);
        indexNode rootKmeans = new indexNode(itemMatrix[0].length);
        setWeight(itemMatrix[0].length, null);// set all as 1
        int depth = 0;
        if (rows.length > leafThreshold && depth < maxDepth) {
            createChildren(root, leafThreshold, depth + 1, maxDepth);
        }
        root.traverseConvert(rootKmeans, itemMatrix[0].length);
        return rootKmeans;
    }
    
    public static indexNode create(double[][] itemMatrix, int leafThreshold, int maxDepth, double []weightinput, int dimension) {
        int[] rows = new int[itemMatrix.length];
        for (int row = 0; row < itemMatrix.length; row++) {
            rows[row] = row;
        }
        Ball root = new Ball(rows, itemMatrix);
        indexNode rootKmeans = new indexNode(dimension);
        setWeight(dimension, weightinput);
        int depth = 0;
<<<<<<< HEAD
//        split conditions to test the function of each condition
//        if (rows.length > leafThreshold && depth < maxDepth) {
//            createChildren(root, leafThreshold, depth + 1, maxDepth);
//        }
        if (rows.length > leafThreshold){
=======
        if (rows.length > leafThreshold && depth < maxDepth) {
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
            createChildren(root, leafThreshold, depth + 1, maxDepth);
        }
        root.traverseConvert(rootKmeans, dimension);
        return rootKmeans;
    }

    private static void createChildren(Ball parent, int leafThreshold, int depth, int maxDepth) {
        IntArrayList leftRows = new IntArrayList();
        IntArrayList rightRows = new IntArrayList();

        splitItems(parent.getRows(), parent.getItemMatrix(), leftRows, rightRows);
        parent.clearRows();

        Ball leftChild = new Ball(leftRows.toIntArray(), parent.getItemMatrix());
        parent.setLeftChild(leftChild);
<<<<<<< HEAD
//        split conditions to test the function of each condition
//        if (leftChild.getRows().length > leafThreshold && depth < maxDepth) {
//            createChildren(leftChild, leafThreshold, depth + 1, maxDepth);
//        }
        if (leftChild.getRows().length > leafThreshold){
=======
        if (leftChild.getRows().length > leafThreshold && depth < maxDepth) {
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
            createChildren(leftChild, leafThreshold, depth + 1, maxDepth);
        }

        Ball rightChild = new Ball(rightRows.toIntArray(), parent.getItemMatrix());
        parent.setRightChild(rightChild);
<<<<<<< HEAD
//        if (rightChild.getRows().length > leafThreshold) {
//            createChildren(rightChild, leafThreshold, depth + 1, maxDepth);
//        }
        if (rightChild.getRows().length > leafThreshold){
=======
        if (rightChild.getRows().length > leafThreshold) {
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
            createChildren(rightChild, leafThreshold, depth + 1, maxDepth);
        }
    }

    protected static void splitItems(int[] rows, double[][] itemMatrix, IntArrayList leftRows, IntArrayList rightRows) {
        // pick random element
        double[] x = itemMatrix[rows[random.nextInt(rows.length)]];
        // select furthest point A to x
        double[] A = x;
        double dist1 = 0;
        for (int row : rows) {
            double[] y = itemMatrix[row];
            double dist2 = distance2(x, y);
            if (dist2 > dist1) {
                A = y;
                dist1 = dist2;
            }
        }
        // select furthest point B to A
        double[] B = A;
        dist1 = 0;
        for (int row : rows) {
            double[] y = itemMatrix[row];
            double dist2 = distance2(A, y);
            if (dist2 > dist1) {
                B = y;
                dist1 = dist2;
            }
        }

        // split data according to A and B proximity
        for (int row : rows) {
            double[] y = itemMatrix[row];
            double distA = distance2(A, y);
            double distB = distance2(B, y);

            if (distA <= distB) {
                leftRows.add(row);
            } else {
                rightRows.add(row);
            }
        }
    }
    
    //the above function can be optimized to build a balanced kd-tree, by changing the spliting rules to Median of medians or midvalue
    

    public static class Ball extends NodeBall {

        private double[] center;
        private double radius;
        private int[] rows;//it
        private final double[][] itemMatrix;

        public Ball(int[] rows, double[][] itemMatrix) {
            this.rows = rows;
            this.itemMatrix = itemMatrix;
            calculateCenter();
<<<<<<< HEAD
//            try new center method
//            calculateCenterByMBR();
=======
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
            calculateRadius();
        }

        @Override
        public Ball getParent() {
            return (Ball) super.getParent();
        }

        @Override
        public Ball getLeftChild() {
            return (Ball) super.getLeftChild();
        }

        @Override
        public Ball getRightChild() {
            return (Ball) super.getRightChild();
        }

        private void calculateCenter() {
            center = new double[itemMatrix[0].length];
<<<<<<< HEAD
//            double[][] temp = new double[rows.length][2];
            /*if(rows.length == 17) {
                for (int i = 0; i < rows.length; i++) {
                    if (rows[i] == 261) {
                        System.out.println("wtf?");
                    }
                    temp[i] = itemMatrix[rows[i]];
                }
            }*/
=======

>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
            for (int row : rows) {
                for (int i = 0; i < center.length; i++) {
                    center[i] += itemMatrix[row][i];
                }
            }
            for (int i = 0; i < center.length; i++) {
                center[i] /= rows.length;
            }
        }

<<<<<<< HEAD
//        a new method to calculate pivot/center using mbr
        private void calculateCenterByMBR() {
            int dim = itemMatrix[0].length;
            center = new double[dim];
            for (int i = 0; i < dim; i++) {
                double max = Double.MIN_VALUE;
                double min = Double.MAX_VALUE;
                for (int row : rows) {
                    max = Math.max(max, itemMatrix[row][i]);
                    min = Math.min(min, itemMatrix[row][i]);
                }
                center[i] = min + (max - min)/2;
            }
        }

=======
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
        private void calculateRadius() {
            radius = Double.NEGATIVE_INFINITY;

            for (int row : rows) {
                radius = max(radius, distance2(center, itemMatrix[row]));
            }
<<<<<<< HEAD
//            radius = sqrt(radius);
=======
            radius = sqrt(radius);
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
        }

        public double mip(double[] q) {
            return dotProduct(q, center) + radius * norm(q);
        }

        public double mip(Ball ball) {
            double[] p0 = center;
            double[] q0 = ball.getCenter();
            double rp = radius;
            double rq = ball.getRadius();
            return dotProduct(p0, q0) + rp * rq + rq * norm(p0) + rp * norm(q0);
        }

        public double[] getCenter() {
            return center;
        }

        public double getRadius() {
            return radius;
        }

        public int[] getRows() {
            return rows;
        }

        public void clearRows() {
            rows = null;
        }

        public double[][] getItemMatrix() {
            return itemMatrix;
        }
        
        public int traverseConvert(indexNode rootKmeans, int dimension) {
    		rootKmeans.setRadius(radius);    		
    		rootKmeans.setPivot(center);//  		
    		if(rows != null){//for the leaf node
    			Set<Integer> aIntegers = new HashSet<Integer>();
    			double []sumOfPoints = new double[dimension];
    			for(int id : rows) {
<<<<<<< HEAD
//    				aIntegers.add(id+1);// the pointid
                    aIntegers.add(id);
//                    test
                    /*if (aIntegers.contains(itemMatrix.length)) {
                        System.out.println("WTF?");
                    }*/
    				for(int i=0; i<dimension; i++)
    					sumOfPoints[i] += itemMatrix[id][i];
    			}
                /*if (rows.length == 17) {
                    for (int r : rows) {
                        if (r == 261) {
                            System.out.println("RUA!!!");
                        }
                    }
                }*/
=======
    				aIntegers.add(id+1);// the pointid
    				for(int i=0; i<dimension; i++)
    					sumOfPoints[i] += itemMatrix[id][i];
    			}
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
    			rootKmeans.setSum(sumOfPoints);
    			rootKmeans.addPoint(aIntegers);		
    			rootKmeans.setTotalCoveredPoints(aIntegers.size());
    			return aIntegers.size();
    		}else {   			
    			int count = 0;
				indexNode childleftnodekmeans = new indexNode(dimension);
				rootKmeans.addNodes(childleftnodekmeans);
				count += getLeftChild().traverseConvert(childleftnodekmeans, dimension);
				
				indexNode childrightnodekmeans = new indexNode(dimension);
				rootKmeans.addNodes(childrightnodekmeans);
				count += getRightChild().traverseConvert(childrightnodekmeans, dimension);
    			rootKmeans.setTotalCoveredPoints(count);
    			return count;
    		}		
    	}
    }

    public static double distance(double[] x, double[] y) {
        return sqrt(distance2(x, y));
    }

    public static double distance2(double[] x, double[] y) {
        double d = 0.0;
        for (int i = 0; i < x.length; i++) {
        	if(weight==null)//normal case
        		d += (x[i] - y[i]) * (x[i] - y[i]);
        	else
        		d += (x[i] - y[i]) * (x[i] - y[i])*weight[i]*weight[i];
        }
<<<<<<< HEAD
        return sqrt(d);
=======
        return d;
>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
    }

    public static double norm(double[] x) {
        return sqrt(norm2(x));
    }

    public static double norm2(double[] x) {
        return dotProduct(x, x);
    }

    public static double dotProduct(double[] x, double[] y) {
        double p = 0.0;
        for (int i = 0; i < x.length; i++) {
            p += x[i] * y[i];
        }
        return p;
    }
}
