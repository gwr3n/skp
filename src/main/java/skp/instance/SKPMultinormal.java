package skp.instance;

import java.util.Arrays;

import umontreal.ssj.probdistmulti.MultiNormalDist;

public class SKPMultinormal extends SKP {

   double[][] covarianceWeights;
   double capacity;
   
   public SKPMultinormal(double[] expectedValues, double[] expectedWeights, double[][] covarianceWeights, double capacity, double shortageCost) {
      super(expectedValues, expectedWeights, shortageCost);
      this.covarianceWeights = covarianceWeights;
      this.capacity = capacity;
      
      generateInstanceID();
   }
   
   public SKPMultinormal(double[] expectedValues, MultiNormalDist weights, double capacity, double shortageCost) {
      this(expectedValues, weights.getMu(), weights.getCovariance(), capacity, shortageCost);
   }
   
   public SKPMultinormal(double[] expectedValues, double[] expectedWeights, double coefficientOfVariation, double correlationCoefficient, double capacity, double shortageCost) {
      this(expectedValues, expectedWeights, calculateCovariance(expectedWeights, coefficientOfVariation, correlationCoefficient), capacity, shortageCost);
   }
   
   public SKPMultinormal(double[] expectedValues, double[] expectedWeights, double coefficientOfVariation, double capacity, double shortageCost) {
      this(expectedValues, expectedWeights, calculateCovariance(expectedWeights, coefficientOfVariation, 0), capacity, shortageCost);
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
    * Compute covariance according to the special structure $\rho^{|j-i|}\sigma_i\sigma_j$ discussed in [1],
    * which ensures P(d_t=x|d_{t-1}=y) = P(d_t=x|d_{t-1}=y,d_{t-2}=z,...)
    * 
    * [1] M. Xiang, R. Rossi, B. Martin-Barragan, S. A. Tarim, "<a href="https://doi.org/10.1016/j.ejor.2022.04.011">
    * A mathematical programming-based solution method for the nonstationary inventory problem under correlated demand</a>,
    * " European Journal of Operational Research, Elsevier, Vol. 304(2): 515–524, 2023
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
      return Arrays.hashCode(new int[] {Arrays.hashCode(this.expectedValues),
                                        Arrays.hashCode(this.expectedWeights),
                                        Arrays.deepHashCode(this.covarianceWeights),
                                        Double.hashCode(this.capacity), 
                                        Double.hashCode(this.shortageCost)});
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPMultinormal) {
         SKPMultinormal o = (SKPMultinormal) obj;
         return Arrays.equals(this.expectedValues, o.expectedValues) &&
                Arrays.equals(this.expectedWeights, o.expectedWeights) &&
                Arrays.deepEquals(this.covarianceWeights, o.covarianceWeights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
   
   public static SKPMultinormal getTestInstance() {
      double[] expectedValues = {111, 111, 21, 117, 123, 34, 3, 121, 112, 12};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      double rho = 0.5;
      int capacity = 100;
      int shortageCost = 10;
      return new SKPMultinormal(expectedValues, expectedWeights, cv, rho, capacity, shortageCost);
   }
   
   public static SKPMultinormal getTestInstanceUncorrelated() {
      double[] expectedValues = {111, 111, 21, 117, 123, 34, 3, 121, 112, 12};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      double rho = 0;
      int capacity = 100;
      int shortageCost = 10;
      return new SKPMultinormal(expectedValues, expectedWeights, cv, rho, capacity, shortageCost);
   }
   
   /**
    * Generate a test instance whose covariance matrix takes the special structure 
    * $\rho^{|j-i|}\sigma_i\sigma_j$ discussed in [1], which ensures 
    * $P(d_t=x|d_{t-1}=y) = P(d_t=x|d_{t-1}=y,d_{t-2}=z,...)$.<br>
    * <br>
    * [1] M. Xiang, R. Rossi, B. Martin-Barragan, S. A. Tarim, "<a href="https://doi.org/10.1016/j.ejor.2022.04.011">
    * A mathematical programming-based solution method for the nonstationary inventory problem under correlated demand</a>,
    * " European Journal of Operational Research, Elsevier, Vol. 304(2): 515–524, 2023
    */
   public static SKPMultinormal getTestInstanceSpecialStructure() {
      double[] expectedValues = {111, 111, 21, 117, 123, 34, 3, 121, 112, 12};
      double[] expectedWeights = {7,5,9,8,6,5,5,9,7,8};
      double cv = 0.2;
      double rho = 0.9;
      int capacity = 30;
      int shortageCost = 10;
      return new SKPMultinormal(expectedValues, expectedWeights, calculateCovarianceSpecialStructure(expectedWeights, cv, rho), capacity, shortageCost);
   }
}
