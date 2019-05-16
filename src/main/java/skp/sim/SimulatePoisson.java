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
import skp.instance.SKPPoisson;
import skp.milp.SKPPoissonMILP;
import skp.milp.instance.SKPPoissonMILPSolvedInstance;
import skp.utililities.gson.GSONUtility;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.randvar.UniformGen;

public class SimulatePoisson extends Simulate {
   
   SKPPoisson instance;
   
   public SimulatePoisson(SKPPoisson instance) {
      this.instance = instance;
   }
   
   public double simulate(int[] knapsack, int nbSamples) {
      double knapsackValue = 0;
      for(int i = 0; i < knapsack.length; i++) {
         if(knapsack[i] == 1) knapsackValue += this.instance.getExpectedValuesPerUnit()[i]*this.instance.getWeights()[i].getMean(); 
      }
      double[][] sampleMatrix = sampleWeights(knapsack, nbSamples);
      knapsackValue -= Arrays.stream(sampleMatrix)
                             .mapToDouble(row -> instance.getShortageCost()*Math.max(0, Arrays.stream(row)
                                                                                              .sum() - instance.getCapacity()))
                             .sum()/nbSamples;
      return knapsackValue;
   }
   
   private double[][] sampleWeights(int[] knapsack, int nbSamples){
      PoissonDist[] weights = this.instance.getWeights();
      PoissonDist[] reducedWeights = IntStream.iterate(0, i -> i + 1)
                                             .limit(weights.length)
                                             .filter(i -> knapsack[i] == 1)
                                             .mapToObj(i -> weights[i])
                                             .toArray(PoissonDist[]::new);
      
      this.randGenerator.resetStartStream();
      double[][] sampleMatrix = new double[nbSamples][reducedWeights.length];
      for(int i = 0; i < sampleMatrix.length; i++){
         for(int j = 0; j < sampleMatrix[i].length; j++){
            sampleMatrix[i][j] = reducedWeights[j].inverseF(UniformGen.nextDouble(this.randGenerator, 0, 1));
         }
      }
      return sampleMatrix;
   }
   
   public static void main(String args[]) {
      
      SKPPoisson instance = SKPPoisson.getTestInstance();
      
      int partitions = 10;
      int linearizationSamples = 50000;
      int simulationSamples = 100000;
      
      try {
         SKPPoissonMILP milp = new SKPPoissonMILP(instance, partitions, linearizationSamples);
         SKPPoissonMILPSolvedInstance solved = milp.solve(simulationSamples);
         System.out.println(GSONUtility.<SKPPoissonMILPSolvedInstance>printInstanceAsJSON(solved));
      } catch (IloException e) {
         e.printStackTrace();
      }
   }
}
