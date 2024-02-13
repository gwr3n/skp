package skp.sim;

import java.util.Arrays;
import java.util.stream.IntStream;

import skp.instance.SKPGenericDistribution;
import skp.milp.SKPGenericDistributionMILP;
import skp.milp.instance.SKPGenericDistributionMILPSolvedInstance;
import skp.utililities.gson.GSONUtility;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.randvar.UniformGen;

public class SimulateGenericDistribution extends Simulate {

   SKPGenericDistribution instance;
   
   public SimulateGenericDistribution(SKPGenericDistribution instance) {
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
      Distribution[] weights = this.instance.getWeights();
      Distribution[] reducedWeights = IntStream.iterate(0, i -> i + 1)
                                               .limit(weights.length)
                                               .filter(i -> knapsack[i] == 1)
                                               .mapToObj(i -> weights[i])
                                               .toArray(Distribution[]::new);
      
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
      
      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstance();
      
      int linearizationSamples = 100000;
      int simulationSamples = 100000;
      
      try {
         SKPGenericDistributionMILP milp = new SKPGenericDistributionMILP(instance, linearizationSamples);
         SKPGenericDistributionMILPSolvedInstance solved = milp.solve(simulationSamples);
         System.out.println(GSONUtility.<SKPGenericDistributionMILPSolvedInstance>printInstanceAsJSON(solved));
      } catch (Exception e) {
         e.printStackTrace();
      }
   }
   
}
