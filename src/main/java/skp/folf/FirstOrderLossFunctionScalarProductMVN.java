package skp.folf;

import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdistmulti.MultiNormalDist;

public class FirstOrderLossFunctionScalarProductMVN {
   MultiNormalDist distribution;
   
   /**
    * Return the FOLF of the scalar product between the Multivariate Normal distribution and vector X
    * 
    * @param distributions
    * @param X
    */
   public FirstOrderLossFunctionScalarProductMVN(MultiNormalDist distribution){
      this.distribution = distribution;
   }
   
   public double getComplementaryFirstOrderLossFunctionValue(double y, double[] x){
      double mu = 0;
      for(int i = 0; i < distribution.getMean().length; i++) {
         mu += distribution.getMu(i)*x[i];
      }
      
      double variance = 0;
      for(int j = 0; j < distribution.getMean().length; j++) {
         double acc = 0;
         for(int i = 0; i < distribution.getMean().length; i++) {
            acc += x[i]*distribution.getCovariance()[i][j];
         }
         acc *= x[j];
         variance += acc;
      }
      double sigma = Math.sqrt(variance);
      
      return sigma * NormalDist.density01((y - mu)/sigma) + (y - mu) * NormalDist.cdf01((y - mu)/sigma);
   }
   
   public double getFirstOrderLossFunctionValue(double y, double[] x){
      double mu = 0;
      for(int i = 0; i < distribution.getMean().length; i++) {
         mu += distribution.getMu(i)*x[i];
      }
      return this.getComplementaryFirstOrderLossFunctionValue(y, x) - (y - mu);
   }
}
