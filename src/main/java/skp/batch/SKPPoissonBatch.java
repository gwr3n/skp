package skp.batch;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.stream.IntStream;

import com.google.gson.Gson;
import com.google.gson.JsonIOException;
import com.google.gson.stream.JsonReader;

import skp.SKPPoisson;
import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.randvar.UniformIntGen;
import umontreal.ssj.rng.MRG32k3aL;

public class SKPPoissonBatch {
   long[] seed = {1,2,3,4,5,6};
   private MRG32k3aL randGenerator;
   
   public SKPPoissonBatch(){
      this.randGenerator = new MRG32k3aL();
      this.randGenerator.setSeed(seed);
   }
   
   public static void main(String args[]) {
      testMultipleInstances();
   }
   
   /*
    * Multiple instances
    */
   
   public static void testMultipleInstances() {
      SKPPoissonBatch batch = new SKPPoissonBatch();
      int numberOfInstances = 10;
      int instanceSize = 5;
      SKPPoisson[] instances = batch.generateInstances(numberOfInstances, instanceSize);
      String fileName = "scrap/instances.json";
      saveInstanceArrayToGSON(instances, fileName);
      SKPPoisson[] newInstances = retrieveInstanceArray(fileName);
      assert(Arrays.deepEquals(instances, newInstances));
   }
   
   public SKPPoisson[] generateInstances(int numberOfInstances, int instanceSize){
      SKPPoisson[] instances = IntStream.iterate(0, i -> i + 1)
                                        .limit(numberOfInstances)
                                        .mapToObj(i -> new SKPPoisson(
                                              (new RandomVariateGen(this.randGenerator, new UniformDist(15,70))).nextArrayOfDouble(instanceSize),
                                              (new RandomVariateGen(this.randGenerator, new UniformDist(15,70))).nextArrayOfDouble(instanceSize),
                                              UniformIntGen.nextInt(this.randGenerator, 50, 150),
                                              UniformGen.nextDouble(this.randGenerator, 50, 150)))
                                        .toArray(SKPPoisson[]::new);
      return instances;
   }
   
   static void saveInstanceArrayToGSON(SKPPoisson[] instance, String fileName){
      Gson gson = new Gson();
      try {
         FileWriter fw = new FileWriter(fileName);
         gson.toJson(instance, fw);
         fw.close();
      } catch (JsonIOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      } catch (IOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
   
   static SKPPoisson[] retrieveInstanceArray(String fileName){
      Gson gson = new Gson();
      JsonReader reader;
      try {
         FileReader fr = new FileReader(fileName);
         reader = new JsonReader(fr);
         SKPPoisson[] instances = gson.fromJson(reader, SKPPoisson[].class);
         fr.close();
         return instances;
      } catch (IOException e) {
         e.printStackTrace();
         return null;
      }
   }
   
   /*
    * Single instance
    */
   
   static void testSingleInstance() {
      SKPPoisson instance = generateInstance();
      String fileName = "scrap/instance.json";
      SKPPoissonBatch.saveInstanceToGSON(instance, fileName);
      SKPPoisson newInstance = retrieveInstance(fileName);
      assert(newInstance.equals(instance));
   }
   
   static SKPPoisson generateInstance(){
      double[] expectedValuesPerUnit = {2.522727273, 2.642857143, 0.287671233, 7.8, 1.732394366, 2.833333333, 0.230769231, 8.642857143, 4.869565217, 0.8};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      int capacity = 100;
      int shortageCost = 100;
      return new SKPPoisson(expectedValuesPerUnit, expectedWeights, capacity, shortageCost);
   }
   
   static void saveInstanceToGSON(SKPPoisson instance, String fileName){
      Gson gson = new Gson();
      try {
         FileWriter fw = new FileWriter(fileName);
         gson.toJson(instance, fw);
         fw.close();
      } catch (JsonIOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      } catch (IOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
   
   static SKPPoisson retrieveInstance(String fileName){
      Gson gson = new Gson();
      JsonReader reader;
      try {
         FileReader fr = new FileReader(fileName);
         reader = new JsonReader(fr);
         SKPPoisson instance = gson.fromJson(reader, SKPPoisson.class);
         fr.close();
         return instance;
      } catch (IOException e) {
         e.printStackTrace();
         return null;
      }
   }
}
