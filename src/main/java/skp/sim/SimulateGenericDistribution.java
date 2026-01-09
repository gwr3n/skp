package skp.sim;

import java.util.Arrays;
import java.util.stream.IntStream;

import skp.instance.SKPGenericDistribution;
import skp.milp.SKPGenericDistributionLazyCuts;
import skp.milp.instance.SKPGenericDistributionCutsSolvedInstance;
import skp.utilities.gson.GSONUtility;
import skp.utilities.probability.SampleFactory;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.rng.RandomStream;

public class SimulateGenericDistribution extends Simulate {

   SKPGenericDistribution instance;
   
   public SimulateGenericDistribution(SKPGenericDistribution instance) {
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
                             .mapToDouble(row -> instance.getShortageCost()*Math.max(0, Arrays.stream(row)
                                                                                              .sum() - instance.getCapacity()))
                             .sum()/nbSamples;
      return knapsackValue;
   }
   
   public double[] simulateMeanVariance(int[] knapsack, int nbSamples, RandomStream randGenerator) {
      double knapsackValue = 0;
      for(int i = 0; i < knapsack.length; i++) {
         if(knapsack[i] == 1) knapsackValue += this.instance.getExpectedValues()[i]; 
      }
      double[][] sampleMatrix = sampleWeights(knapsack, nbSamples, randGenerator);
      double[] knapsackValues = Arrays.stream(sampleMatrix)
                                      .mapToDouble(row -> instance.getShortageCost()*Math.max(0, Arrays.stream(row)
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
   
   private double[][] sampleWeights(int[] knapsack, int nbSamples, RandomStream rng){
      Distribution[] weights = this.instance.getWeights();
      Distribution[] reducedWeights = IntStream.iterate(0, i -> i + 1)
                                               .limit(weights.length)
                                               .filter(i -> knapsack[i] == 1)
                                               .mapToObj(i -> weights[i])
                                               .toArray(Distribution[]::new);
      
      rng.resetNextSubstream();
      double[][] sampleMatrix;
      switch(Simulate.samplingStrategy) {
      case LHS:
         sampleMatrix = SampleFactory.getNextLHSample(reducedWeights, nbSamples, rng);
         break;
      case SRS:
      default:
         sampleMatrix = SampleFactory.getNextSimpleRandomSample(reducedWeights, nbSamples, rng);
      }
      return sampleMatrix;
   }
   
   public static void main(String args[]) {
      
      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstanceLarge();
      
      int linearizationSamples = 10000;
      int simulationRuns = 100000;
      
      try {
         SKPGenericDistributionLazyCuts milp = new SKPGenericDistributionLazyCuts(instance, linearizationSamples, simulationRuns);
         SKPGenericDistributionCutsSolvedInstance solved = milp.solve();
         System.out.println(GSONUtility.<SKPGenericDistributionCutsSolvedInstance>printInstanceAsJSON(solved));
      } catch (Exception e) {
         e.printStackTrace();
      }
   }
   
}
