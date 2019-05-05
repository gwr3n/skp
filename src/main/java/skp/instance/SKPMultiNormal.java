package skp.instance;

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
   
   public SKPMultiNormal(double[] expectedValuesPerUnit, double[] expectedWeights, double coefficientOfVariation, double correlationCoefficient, double capacity, double shortageCost) {
      this(expectedValuesPerUnit, expectedWeights, calculateCovariance(expectedWeights, coefficientOfVariation, correlationCoefficient), capacity, shortageCost);
   }
   
   public static double[][] calculateCovariance(double [] means, double cv, double rho){
      double[] stdDemand =new double [means.length];
      for (int i = 0; i < means.length; i ++) {
         stdDemand[i] = cv * means[i];
      }
      
      double[][] covariance = new double [means.length][means.length];
      for (int row=0; row<covariance.length;row++) {
         for (int col=0; col<covariance[row].length;col++) {
            if (row==col) {
               covariance[row][col]=stdDemand[row]*stdDemand[col];
            } else {
               covariance[row][col]=stdDemand[row]*stdDemand[col]*rho;
            }
         }
      }
      return covariance;
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
   
   public static SKPMultiNormal getTestInstance() {
      double[] expectedValuesPerUnit = {2.522727273, 2.642857143, 0.287671233, 7.8, 1.732394366, 2.833333333, 0.230769231, 8.642857143, 4.869565217, 0.8};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      double rho = 0.5;
      int capacity = 100;
      int shortageCost = 100;
      return new SKPMultiNormal(expectedValuesPerUnit, expectedWeights, cv, rho, capacity, shortageCost);
   }
}
