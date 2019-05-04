package skp;

import java.math.BigInteger;
import java.util.Arrays;

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
      this.expectedValuesPerUnit = expectedValuesPerUnit;
      this.expectedWeights = weights.getMu();
      this.covarianceWeights = weights.getCovariance();
      this.capacity = capacity;
      this.shortageCost = shortageCost;
      
      generateInstanceID();
   }
   
   private void generateInstanceID() {
      String intHash = ""+this.hashCode();
      BigInteger hashCode = new BigInteger(intHash);
      instanceID = hashCode.toString(16);
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
