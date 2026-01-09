package skp.utilities.probability;

import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.apache.commons.math3.linear.MatrixUtils;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.randvar.UniformIntGen;
import umontreal.ssj.rng.RandomStream;

public class SampleFactory {
   
   /**
    * Implements Simple Random Sampling
    * @param distributions array of distributions to be sampled
    * @return a Simple Random Sample for the distributions in {@code distributions}
    */
   public static double[][] getNextSimpleRandomSample(Distribution[] distributions, int samples, RandomStream stream){
      double[][] sampleMatrix = new double[samples][distributions.length];
      for(int i = 0; i < sampleMatrix.length; i++){
         for(int j = 0; j < sampleMatrix[i].length; j++){
            sampleMatrix[i][j] = distributions[j].inverseF(UniformGen.nextDouble(stream, 0, 1));
         }
      }
      return sampleMatrix;
   }
   
   /**
    * Implements Latin Hypercube Sampling as originally introduced in 
    * 
    * McKay, M.D.; Beckman, R.J.; Conover, W.J. (May 1979). 
    * "A Comparison of Three Methods for Selecting Values of Input Variables 
    * in the Analysis of Output from a Computer Code". 
    * Technometrics 21 (2): 239–245.
    * 
    * @param distributions array of distributions to be sampled 
    * @param samples number of samples
    * @return a Latin Hypercube Sample for the distributions in {@code distributions}
    */
   public static double[][] getNextLHSample(Distribution[] distributions, int samples, RandomStream stream){
      if(distributions.length == 0)
         return new double[samples][distributions.length];
      double x[][] = new double[distributions.length][samples];
      x = IntStream.iterate(0, d -> d + 1)
                   .limit(distributions.length)
                   .mapToObj(
                         d -> DoubleStream.iterate(0, i -> i + 1.0/samples)
                                          .limit(samples)
                                          .map(i -> distributions[d].inverseF(i + UniformGen.nextDouble(stream, 0, 1.0/samples)))
                                          .toArray())
                   .toArray(double[][]::new);   
      for(int i = 0; i < x.length; i++){
         shuffle(x[i], stream);
      }
      double[][] xT = MatrixUtils.createRealMatrix(x).transpose().getData();
      return xT;
   }
   
   /**
    * Returns a random shuffle of {@code sample}.
    * 
    * @param sample the original sample.
    * @return a random shuffle of {@code sample}.
    */
   private static double[] shuffle(double[] sample, RandomStream stream){
      for(int i = 0; i < sample.length; i++){
         int j = UniformIntGen.nextInt(stream, 0, sample.length - 1);
         double temp = sample[i];
         sample[i] = sample[j];
         sample[j] = temp;
      }
      return sample;
   }
}
