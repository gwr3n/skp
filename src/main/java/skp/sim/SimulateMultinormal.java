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

import skp.instance.SKPMultinormal;
import skp.milp.PWAPPROXIMATION;
import skp.milp.SKPMultinormalMILP;
import skp.milp.instance.SKPMultinormalMILPSolvedInstance;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalGen;
import umontreal.ssj.randvarmulti.MultinormalPCAGen;
import umontreal.ssj.rng.RandomStream;

public class SimulateMultinormal extends Simulate {
   
   SKPMultinormal instance;
   
   public SimulateMultinormal(SKPMultinormal instance) {
      super();
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
   
   public double[] simulateMeanVariance(int[] knapsack, int nbSamples) {
      double knapsackValue = 0;
      for(int i = 0; i < knapsack.length; i++) {
         if(knapsack[i] == 1) knapsackValue += this.instance.getExpectedValues()[i]; 
      }
      double[][] sampleMatrix = sampleWeights(nbSamples);
      double[] knapsackValues = Arrays.stream(sampleMatrix)
                                      .mapToDouble(row -> instance.getShortageCost()*Math.max(0, IntStream.iterate(0, i -> i + 1)
                                                                                                          .limit(row.length)
                                                                                                          .filter(i -> knapsack[i] == 1)
                                                                                                          .mapToDouble(i -> row[i])
                                                                                                          .sum() - instance.getCapacity()))
                                      .toArray();
      for(int s = 0; s < knapsackValues.length; s++)
         knapsackValues[s] = knapsackValue - knapsackValues[s];
      
      return new double[] {calculateMean(knapsackValues), calculateVariance(knapsackValues)};
   }
   
   public static double calculateVariance(double[] numbers) {
      double mean = calculateMean(numbers);
      double sum = 0;
      for (double num : numbers) {
          sum += Math.pow(num - mean, 2);
      }
      return sum / numbers.length;
   }

   public static double calculateMean(double[] numbers) {
      double sum = 0;
      for (double num : numbers) {
          sum += num;
      }
      return sum / numbers.length;
   }
   
   private double[][] sampleWeights(int nbSamples) {
      double[] mu = this.instance.getWeights().getMean();
      double[][] sigma = this.instance.getWeights().getCovariance();
      
      this.randGenerator.resetStartStream();
      NormalGen standardNormal = new NormalGen(this.randGenerator, 0, 1);
      //MultinormalGen gen = new MultinormalCholeskyGen(standardNormal, mu, sigma);
      MultinormalGen gen = new MultinormalPCAGen(standardNormal, mu, sigma);
      double[][] points = new double[nbSamples][this.instance.getItems()];
      for(int i = 0; i < nbSamples; i++)
         gen.nextPoint(points[i]);
      return points;
   }
   
   public double[][] sampleWeights(int nbSamples, RandomStream rng) {
      double[] mu = this.instance.getWeights().getMean();
      double[][] sigma = this.instance.getWeights().getCovariance();
      
      rng.resetStartStream();
      NormalGen standardNormal = new NormalGen(rng, 0, 1);
      //MultinormalGen gen = new MultinormalCholeskyGen(standardNormal, mu, sigma);
      MultinormalGen gen = new MultinormalPCAGen(standardNormal, mu, sigma);
      double[][] points = new double[nbSamples][this.instance.getItems()];
      for(int i = 0; i < nbSamples; i++)
         gen.nextPoint(points[i]);
      return points;
   }
   
   public static void main(String args[]) {
      
      SKPMultinormal instance = SKPMultinormal.getTestInstanceSpecialStructure();
      
      int partitions = 10;
      int simulationRuns = 1000000;
      
      try {
         SKPMultinormalMILP milp = new SKPMultinormalMILP(instance, partitions, PWAPPROXIMATION.EDMUNDSON_MADANSKI);
         
         System.out.println(GSONUtility.<SKPMultinormalMILPSolvedInstance>printInstanceAsJSON(milp.solve(simulationRuns)));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
   }

}
