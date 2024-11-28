package skp.folf;

import java.util.Arrays;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.EmpiricalDist;

public class FirstOrderLossFunctionScalarProduct extends FirstOrderLossFunction{
   /**
    * Return the FOLF of the scalar product between the variables in distributions and vector X
    * 
    * @param distributions
    * @param X
    */
   public FirstOrderLossFunctionScalarProduct(Distribution[] distributions){
      super(distributions);
   }
   
   public EmpiricalDist getEmpiricalDistribution(int nbSamples, double[] x){
      double[][] sampleMatrix = this.sample(nbSamples);
      double[] observations = new double[nbSamples];
      for(int i = 0; i < sampleMatrix.length; i++){
         for(int j = 0; j < sampleMatrix[i].length; j++){
            observations[i] += sampleMatrix[i][j]*x[j];
         }
      }
      Arrays.sort(observations);
      EmpiricalDist empDistribution = new EmpiricalDist(observations);
      return empDistribution;
   }
   
   public double getComplementaryFirstOrderLossFunctionValue(double y, double[] x, int nbSamples){
      EmpiricalDist empDistribution = this.getEmpiricalDistribution(nbSamples, x);
      double value = 0;
      for(int i = 0; i < empDistribution.getN(); i++){
         value += Math.max(y-empDistribution.getObs(i),0)/empDistribution.getN();
      }
      return value;
   }
   
   public double getFirstOrderLossFunctionValue(double y, double[] x, int nbSamples){
      EmpiricalDist empDistribution = this.getEmpiricalDistribution(nbSamples, x);
      double value = 0;
      for(int i = 0; i < empDistribution.getN(); i++){
         value += Math.max(empDistribution.getObs(i)-y,0)/empDistribution.getN();
      }
      return value;
   }
}
