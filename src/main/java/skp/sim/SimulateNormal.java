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
import skp.instance.SKPNormal;
import skp.milp.PWAPPROXIMATION;
import skp.milp.SKPNormalMILP;
import skp.milp.instance.SKPNormalMILPSolvedInstance;
import skp.utilities.gson.GSONUtility;
import skp.utilities.probability.SampleFactory;

import umontreal.ssj.probdist.Distribution;

public class SimulateNormal extends Simulate {
   
   SKPNormal instance;
   
   public SimulateNormal(SKPNormal instance) {
      super();
      this.instance = instance;
   }
   
   public double simulate(int[] knapsack, int nbSamples) {
      double knapsackValue = 0;
      for(int i = 0; i < knapsack.length; i++) {
         if(knapsack[i] == 1) knapsackValue += this.instance.getExpectedValues()[i]; 
      }
      double[][] sampleMatrix = sampleWeights(knapsack, nbSamples);
      knapsackValue -= Arrays.stream(sampleMatrix)
                             .mapToDouble(row -> instance.getShortageCost()*Math.max(0, Arrays.stream(row).sum() - instance.getCapacity()))
                             .sum()/nbSamples;
      return knapsackValue;
   }
   
   private double[][] sampleWeights(int[] knapsack, int nbSamples){
      Distribution[] weights = this.instance.getWeights();
      Distribution[] reducedWeights = IntStream.iterate(0, i -> i + 1)
                                             .limit(weights.length)
                                             .filter(i -> knapsack[i] == 1)
                                             .mapToObj(i -> weights[i])
                                             .toArray(Distribution[]::new);
      
      this.randGenerator.resetStartStream();
      double[][] sampleMatrix;
      switch(Simulate.samplingStrategy) {
      case LHS:
         sampleMatrix = SampleFactory.getNextLHSample(reducedWeights, nbSamples, randGenerator);
         break;
      case SRS:
      default:
         sampleMatrix = SampleFactory.getNextSimpleRandomSample(reducedWeights, nbSamples, randGenerator);
      }
      return sampleMatrix;
   }
   
   public static void main(String args[]) {
      
      SKPNormal instance = SKPNormal.getTestInstance();
      
      int partitions = 10;
      int simulationSamples = 100000;
      
      try {
         SKPNormalMILP milp = new SKPNormalMILP(instance, partitions, PWAPPROXIMATION.EDMUNDSON_MADANSKI);
         SKPNormalMILPSolvedInstance solved = milp.solve(simulationSamples);
         System.out.println(GSONUtility.<SKPNormalMILPSolvedInstance>printInstanceAsJSON(solved));
      } catch (IloException e) {
         e.printStackTrace();
      }
   }
}
