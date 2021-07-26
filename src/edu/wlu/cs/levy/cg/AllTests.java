package edu.wlu.cs.levy.cg;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

import edu.wlu.cs.levy.cg.HammingDistanceTest;

@RunWith(Suite.class)
@Suite.SuiteClasses({
	HammingDistanceTest.class,
	KDTreeTest.class,
	NearestNeighborListTest.class
})

public class AllTests {
	//don't have to create or destroy anything
}
