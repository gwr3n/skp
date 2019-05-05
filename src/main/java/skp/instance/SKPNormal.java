package skp.instance;

import java.util.Arrays;
import java.util.stream.IntStream;

import skp.utililities.hash.SHA;
import umontreal.ssj.probdist.NormalDist;

public class SKPNormal {
   String instanceID;
   double[] expectedValuesPerUnit;
   double[] expectedWeights;
   double[] stdWeights;
   double capacity;
   double shortageCost;
   
   public SKPNormal(double[] expectedValuesPerUnit, double[] expectedWeights, double[] stdWeights, double capacity, double shortageCost) {
      this.expectedValuesPerUnit = expectedValuesPerUnit;
      this.expectedWeights = expectedWeights;
      this.stdWeights = stdWeights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
      
      generateInstanceID();
   }
   
   public SKPNormal(double[] expectedValuesPerUnit, double[] expectedWeights, double coefficientOfVariation, double capacity, double shortageCost) {
      this(expectedValuesPerUnit, 
            expectedWeights,
            IntStream.iterate(0, i -> i + 1)
                     .limit(expectedWeights.length)
                     .mapToDouble(i -> expectedWeights[i]*coefficientOfVariation)
                     .toArray(),
            capacity, 
            shortageCost);
   }
   
   public SKPNormal(double[] expectedValuesPerUnit, NormalDist[] weights, double capacity, double shortageCost) {
      this(expectedValuesPerUnit,
            IntStream.iterate(0, i -> i + 1)
                     .limit(weights.length)
                     .mapToDouble(i -> weights[i].getMean())
                     .toArray(),
            IntStream.iterate(0, i -> i + 1)
                     .limit(weights.length)
                     .mapToDouble(i -> weights[i].getSigma())
                     .toArray(),
            capacity,
            shortageCost);
   }
   
   private void generateInstanceID() {
      String intHash = ""+this.hashCode();
      instanceID = SHA.generateSHA256(intHash);
   }
   
   public String getInstanceID() {
      return this.instanceID;
   }
   
   public int getItems() {
      return this.expectedValuesPerUnit.length;
   }
   
   public double getCapacity() {
      return this.capacity;
   }
   
   public double getShortageCost() {
      return this.shortageCost;
   }
   
   public double[] getExpectedValuesPerUnit() {
      return this.expectedValuesPerUnit;
   }
   
   public NormalDist[] getWeights() {
      return IntStream.iterate(0, i -> i + 1).limit(expectedWeights.length)
                      .mapToObj(i -> new NormalDist(expectedWeights[i], stdWeights[i]))
                      .toArray(NormalDist[]::new);
   }
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(this.expectedValuesPerUnit) +
            Arrays.hashCode(this.expectedWeights) +
            Arrays.hashCode(this.stdWeights) +
            Double.hashCode(this.capacity) + 
            Double.hashCode(this.shortageCost);
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPNormal) {
         SKPNormal o = (SKPNormal) obj;
         return Arrays.equals(this.expectedValuesPerUnit, o.expectedValuesPerUnit) &&
                Arrays.equals(this.expectedWeights, o.expectedWeights) &&
                Arrays.equals(this.stdWeights, o.stdWeights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
}
