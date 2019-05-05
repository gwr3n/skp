package skp.instance;

import java.util.Arrays;
import java.util.stream.IntStream;

import skp.utililities.hash.SHA;
import umontreal.ssj.probdist.PoissonDist;

public class SKPPoisson {
   String instanceID;
   double[] expectedValuesPerUnit;
   double[] expectedWeights;
   int capacity;
   double shortageCost;
   
   public SKPPoisson(double[] expectedValuesPerUnit, double[] expectedWeights, int capacity, double shortageCost) {
      this.expectedValuesPerUnit = expectedValuesPerUnit;
      this.expectedWeights = expectedWeights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
      
      generateInstanceID();
   }
   
   private void generateInstanceID() {
      String intHash = ""+this.hashCode();
      instanceID = SHA.generateSHA256(intHash);
   }
   
   public SKPPoisson(double[] expectedValuesPerUnit, PoissonDist[] weights, int capacity, double shortageCost) {
      this.expectedValuesPerUnit = expectedValuesPerUnit;
      this.expectedWeights =  IntStream.iterate(0, i -> i + 1).limit(weights.length)
                                       .mapToDouble(i -> weights[i].getLambda())
                                       .toArray();
      this.capacity = capacity;
      this.shortageCost = shortageCost;
      
      generateInstanceID();
   }
   
   public String getInstanceID() {
      return this.instanceID;
   }
   
   public int getItems() {
      return this.expectedValuesPerUnit.length;
   }
   
   public int getCapacity() {
      return this.capacity;
   }
   
   public double getShortageCost() {
      return this.shortageCost;
   }
   
   public double[] getExpectedValuesPerUnit() {
      return this.expectedValuesPerUnit;
   }
   
   public PoissonDist[] getWeights() {
      return IntStream.iterate(0, i -> i + 1).limit(expectedWeights.length)
                      .mapToObj(i -> new PoissonDist(expectedWeights[i]))
                      .toArray(PoissonDist[]::new);
   }
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(this.expectedValuesPerUnit) +
            Arrays.hashCode(this.expectedWeights) +
            Double.hashCode(this.capacity) + 
            Double.hashCode(this.shortageCost);
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPPoisson) {
         SKPPoisson o = (SKPPoisson) obj;
         return Arrays.equals(this.expectedValuesPerUnit, o.expectedValuesPerUnit) &&
                Arrays.equals(this.expectedWeights, o.expectedWeights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
   
   public static SKPPoisson getTestInstance() {
      double[] expectedValuesPerUnit = {2.522727273, 2.642857143, 0.287671233, 7.8, 1.732394366, 2.833333333, 0.230769231, 8.642857143, 4.869565217, 0.8};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      int capacity = 100;
      int shortageCost = 100;
      return new SKPPoisson(expectedValuesPerUnit, expectedWeights, capacity, shortageCost);
   }
}
