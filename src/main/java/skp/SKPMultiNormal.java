package skp;

import java.util.Arrays;

import skp.utililities.hash.SHA;
import umontreal.ssj.probdistmulti.MultiNormalDist;

public class SKPMultiNormal {
   String instanceID;
   double[] expectedValuesPerUnit;
   double[] expectedWeights;
   double[][] covarianceWeights;
   double capacity;
   double shortageCost;
   
   public SKPMultiNormal(double[] expectedValuesPerUnit, double[] expectedWeights, double[][] covarianceWeights, double capacity, double shortageCost) {
      this.expectedValuesPerUnit = expectedValuesPerUnit;
      this.expectedWeights = expectedWeights;
      this.covarianceWeights = covarianceWeights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
      
      generateInstanceID();
   }
   
   public SKPMultiNormal(double[] expectedValuesPerUnit, MultiNormalDist weights, double capacity, double shortageCost) {
      this(expectedValuesPerUnit, weights.getMu(), weights.getCovariance(), capacity, shortageCost);
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
   
   public MultiNormalDist getWeights() {
      return new MultiNormalDist(this.expectedWeights, this.covarianceWeights);
   }
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(this.expectedValuesPerUnit) +
            Arrays.hashCode(this.expectedWeights) +
            Arrays.deepHashCode(this.covarianceWeights) +
            Double.hashCode(this.capacity) + 
            Double.hashCode(this.shortageCost);
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPMultiNormal) {
         SKPMultiNormal o = (SKPMultiNormal) obj;
         return Arrays.equals(this.expectedValuesPerUnit, o.expectedValuesPerUnit) &&
                Arrays.equals(this.expectedWeights, o.expectedWeights) &&
                Arrays.deepEquals(this.covarianceWeights, o.covarianceWeights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
}
