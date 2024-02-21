package skp.instance;

import java.util.Arrays;

import umontreal.ssj.probdistmulti.MultiNormalDist;

public class SKPMultinormal extends SKP {

   double[][] covarianceWeights;
   double capacity;
   
   public SKPMultinormal(double[] expectedValuesPerUnit, double[] expectedWeights, double[][] covarianceWeights, double capacity, double shortageCost) {
      super(expectedValuesPerUnit, expectedWeights, shortageCost);
      this.covarianceWeights = covarianceWeights;
      this.capacity = capacity;
      
      generateInstanceID();
   }
   
   public SKPMultinormal(double[] expectedValuesPerUnit, MultiNormalDist weights, double capacity, double shortageCost) {
      this(expectedValuesPerUnit, weights.getMu(), weights.getCovariance(), capacity, shortageCost);
   }
   
   public SKPMultinormal(double[] expectedValuesPerUnit, double[] expectedWeights, double coefficientOfVariation, double correlationCoefficient, double capacity, double shortageCost) {
      this(expectedValuesPerUnit, expectedWeights, calculateCovariance(expectedWeights, coefficientOfVariation, correlationCoefficient), capacity, shortageCost);
   }
   
   public SKPMultinormal(double[] expectedValuesPerUnit, double[] expectedWeights, double coefficientOfVariation, double capacity, double shortageCost) {
      this(expectedValuesPerUnit, expectedWeights, calculateCovariance(expectedWeights, coefficientOfVariation, 0), capacity, shortageCost);
   }
   
   public static double[][] calculateCovariance(double [] means, double cv, double rho){
      double[] stdDemand = new double [means.length];
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
   
   /**
    * Compute covariance according to the special structure $\rho^{|j-i|}\sigma_i\sigma_j$
    * which ensures P(d_t=x|d_{t-1}=y) = P(d_t=x|d_{t-1}=y,d_{t-2}=z,...)
    * 
    * see https://doi.org/10.1016/j.ejor.2022.04.011
    */
   public static double[][] calculateCovarianceSpecialStructure(double [] means, double cv, double rho){
      double[] stdDemand = new double [means.length];
      for (int i = 0; i < means.length; i ++) {
         stdDemand[i] = cv * means[i];
      }
      
      double[][] covariance = new double [means.length][means.length];
      for (int row=0; row<covariance.length;row++) {
         for (int col=0; col<covariance[row].length;col++) {
            if (row==col) {
               covariance[row][col]=stdDemand[row]*stdDemand[col];
            } else {
               covariance[row][col]=stdDemand[row]*stdDemand[col]*Math.pow(rho, Math.abs(col-row));
            }
         }
      }
      return covariance;
   }
   
   public double getCapacity() {
      return this.capacity;
   }
   
   public MultiNormalDist getWeights() {
      return new MultiNormalDist(this.expectedWeights, this.covarianceWeights);
   }
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(new int[] {Arrays.hashCode(this.expectedValuesPerUnit),
                                        Arrays.hashCode(this.expectedWeights),
                                        Arrays.deepHashCode(this.covarianceWeights),
                                        Double.hashCode(this.capacity), 
                                        Double.hashCode(this.shortageCost)});
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPMultinormal) {
         SKPMultinormal o = (SKPMultinormal) obj;
         return Arrays.equals(this.expectedValuesPerUnit, o.expectedValuesPerUnit) &&
                Arrays.equals(this.expectedWeights, o.expectedWeights) &&
                Arrays.deepEquals(this.covarianceWeights, o.covarianceWeights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
   
   public static SKPMultinormal getTestInstance() {
      double[] expectedValuesPerUnit = {2.522727273, 2.642857143, 0.287671233, 7.8, 1.732394366, 2.833333333, 0.230769231, 8.642857143, 4.869565217, 0.8};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      double rho = 0.5;
      int capacity = 100;
      int shortageCost = 100;
      return new SKPMultinormal(expectedValuesPerUnit, expectedWeights, cv, rho, capacity, shortageCost);
   }
   
   /**
    * Generate a test instance according to the special structure $\rho^{|j-i|}\sigma_i\sigma_j$
    * 
    * see https://doi.org/10.1016/j.ejor.2022.04.011
    */
   public static SKPMultinormal getTestInstanceSpecialStructure() {
      double[] expectedValuesPerUnit = {2.522727273, 2.642857143, 0.287671233, 7.8, 1.732394366};
      double[] expectedWeights = {44,42,73,15,71};
      double cv = 0.2;
      double rho = 0.9;
      int capacity = 100;
      int shortageCost = 100;
      return new SKPMultinormal(expectedValuesPerUnit, expectedWeights, calculateCovarianceSpecialStructure(expectedWeights, cv, rho), capacity, shortageCost);
   }
}
