/**
 * To run from Mac OS
 * 
 * -Djava.library.path=/Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
 * 
 * Environment variable
 * 
 * DYLD_LIBRARY_PATH /Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
 * 
 * @author gwren
 *
 */

package skp.sim;

import java.util.Arrays;
import java.util.stream.IntStream;

import ilog.concert.IloException;
import skp.SKPMultiNormal;
import skp.milp.SKPMultiNormalMILP;
import umontreal.ssj.probdistmulti.MultiNormalDist;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalCholeskyGen;
import umontreal.ssj.randvarmulti.MultinormalGen;
import umontreal.ssj.rng.MRG32k3aL;

public class SimulateMultiNormal {
   SKPMultiNormal instance;
   private MRG32k3aL randGenerator;
   
   public SimulateMultiNormal(SKPMultiNormal instance, long[] seed) {
      this.randGenerator = new MRG32k3aL();
      this.randGenerator.setSeed(seed);
      this.instance = instance;
   }
   
   public double simulate(int[] knapsack, int nbSamples) {
      double knapsackValue = 0;
      for(int i = 0; i < knapsack.length; i++) {
         if(knapsack[i] == 1) knapsackValue += this.instance.getExpectedValues()[i]; 
      }
      double[][] sampleMatrix = sampleWeights(nbSamples);
      knapsackValue -= Arrays.stream(sampleMatrix)
                             .mapToDouble(row -> instance.getShortageCost()*Math.max(0, IntStream.iterate(0, i -> i + 1)
                                                                                                 .limit(row.length)
                                                                                                 .filter(i -> knapsack[i] == 1)
                                                                                                 .mapToDouble(i -> row[i])
                                                                                                 .sum() - instance.getCapacity()))
                             .sum()/nbSamples;
      return knapsackValue;
   }
   
   public double[][] sampleWeights(int nbSamples) {
      double[] mu = this.instance.getWeights().getMean();
      double[][] sigma = this.instance.getWeights().getCovariance();
      this.randGenerator.resetStartStream();
      NormalGen standardNormal = new NormalGen(this.randGenerator, 0, 1);
      MultinormalGen gen = new MultinormalCholeskyGen(standardNormal, mu, sigma);
      double[][] points = new double[nbSamples][this.instance.getItems()];
      for(int i = 0; i < nbSamples; i++)
         gen.nextPoint(points[i]);
      return points;
   }
   
   public static double[][] calculateCovariance(double [] means, double cv, double rho){
      double [] stdDemand =new double [means.length];
      for (int i = 0; i < means.length; i ++) {
         stdDemand[i] = cv*means[i];
      }
      
      double [][] covariance = new double [means.length][means.length];
      
      for (int row=0; row<covariance.length;row++) {
         for (int col=0; col<covariance[row].length;col++) {
            if (row==col) {
               covariance[row][col]=stdDemand[row]*stdDemand[col];
            } else if (col==row+1 | col==row-1) {
               covariance[row][col]=stdDemand[row]*stdDemand[col]*rho;
            } else  {
               covariance[row][col]=0;
            }
         }
      }
      return covariance;
   }
   
   public static void main(String args[]) {
      long[] seed = {1,2,3,4,5,6};
      
      double[] expectedValues = {111,111,21,117,123,34,3,121,112,12};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      double rho = 0.5;
      double[][] varianceCovarianceWeights = calculateCovariance(expectedWeights, cv, rho);
      
      int capacity = 100;
      int shortageCost = 100;
      
      MultiNormalDist weights = new MultiNormalDist(expectedWeights, varianceCovarianceWeights);
      
      SKPMultiNormal instance = new SKPMultiNormal(expectedValues, weights, capacity, shortageCost);
      
      int partitions = 10;
      SKPMultiNormalMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new SKPMultiNormalMILP(instance, partitions);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      
      SimulateMultiNormal sim = new SimulateMultiNormal(instance, seed);
      double milpSolutionValue = milp.getSolutionValue();
      int nbSamples = 20000;
      double simSolutionValue = sim.simulate(knapsack, nbSamples);
      
      System.out.println("MILP: "+milpSolutionValue);
      System.out.println("Simulation: "+simSolutionValue);
      System.out.println("Linearization gap: "+100*(simSolutionValue-milpSolutionValue)/simSolutionValue);
   }

}
