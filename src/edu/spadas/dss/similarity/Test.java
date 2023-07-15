package edu.nyu.dss.similarity;

import java.util.Arrays;
import java.util.Random;

public class Test {
    public static void main(String[] args) {
        int[][] triplets = generateTriplets(100);

        // 打印生成的三元组
        for (int[] triplet : triplets) {
            System.out.println(Arrays.toString(triplet));
        }
    }

    public static int[][] generateTriplets(int count) {
        int[][] triplets = new int[count][3];

        Random random = new Random();
        for (int i = 0; i < count; i++) {
            triplets[i][0] = random.nextInt(101); // 生成0~100的随机数
            triplets[i][1] = random.nextInt(101);
            triplets[i][2] = random.nextInt(101);
        }

        return triplets;
    }
}
