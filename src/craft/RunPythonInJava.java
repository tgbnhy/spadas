package craft;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class RunPythonInJava {
    public static void main(String[] args) {
        Process proc;
        try {
            proc = Runtime.getRuntime().exec("python \"E:\\temp\\Python Projects\\kneed\\kneed_test.py\"");
            BufferedReader in = new BufferedReader(new InputStreamReader(proc.getInputStream()));
//            String line = null;
//            while ((line = in.readLine()) != null) {
//
//                System.out.println(line);
//            }
            System.out.println(in.readLine());
            System.out.println(in.readLine());
            in.close();
            proc.waitFor();
        } catch (InterruptedException | IOException e) {
            e.printStackTrace();
        }
    }
}
