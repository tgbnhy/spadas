<<<<<<< HEAD
package au.edu.rmit.trajectory.clustering.kpaths;

import java.util.ArrayList;

class Path implements Comparable<Path>{//queue is very slow
    private ArrayList<Integer> pathIDs;
    private double score;
    private double lowerbound;// the sum weight

    public Path(ArrayList<Integer> pathIDs, double Weight, double lowerbound) {
        this.pathIDs = pathIDs;
        this.lowerbound = lowerbound;
        this.score = Weight;
    }
    
    public ArrayList<Integer> getPath() {
    	return pathIDs;
    }
    
    public double getLowerbound() {
    	return lowerbound;
    }
    
    public double getScore() {
    	return score;
    }
    
    @Override
    public int compareTo(Path other) {
        return (int)(this.getLowerbound() - other.getLowerbound());
    }
}

=======
package au.edu.rmit.trajectory.clustering.kpaths;

import java.util.ArrayList;

class Path implements Comparable<Path>{//queue is very slow
    private ArrayList<Integer> pathIDs;
    private double score;
    private double lowerbound;// the sum weight

    public Path(ArrayList<Integer> pathIDs, double Weight, double lowerbound) {
        this.pathIDs = pathIDs;
        this.lowerbound = lowerbound;
        this.score = Weight;
    }
    
    public ArrayList<Integer> getPath() {
    	return pathIDs;
    }
    
    public double getLowerbound() {
    	return lowerbound;
    }
    
    public double getScore() {
    	return score;
    }
    
    @Override
    public int compareTo(Path other) {
        return (int)(this.getLowerbound() - other.getLowerbound());
    }
}

>>>>>>> 76a8efdf061f79d681a7e8054c4dbe70dc82b9d3
