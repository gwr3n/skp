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

import skp.SKPNormal;
import skp.milp.SKPNormalMILP;

import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.rng.MRG32k3aL;

public class SimulateNormal {
   SKPNormal instance;
   private MRG32k3aL randGenerator;
   
   public SimulateNormal(SKPNormal instance, long[] seed) {
      this.randGenerator = new MRG32k3aL();
      this.randGenerator.setSeed(seed);
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
   
   private double[][] sampleWeights(int[] knapsack, int nbSamples){
      NormalDist[] weights = this.instance.getWeights();
      NormalDist[] reducedWeights = IntStream.iterate(0, i -> i + 1)
                                             .limit(weights.length)
                                             .filter(i -> knapsack[i] == 1)
                                             .mapToObj(i -> weights[i])
                                             .toArray(NormalDist[]::new);
      
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
      long[] seed = {1,2,3,4,5,6};
      
      double[] expectedValues = {111,111,21,117,123,34,3,121,112,12};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      double[] varianceWeights = IntStream.iterate(0, i -> i + 1)
                                          .limit(expectedWeights.length)
                                          .mapToDouble(i -> expectedWeights[i]*cv).toArray();
      double[] standardDeviationWeights = Arrays.stream(varianceWeights)
                                                .map(w -> Math.sqrt(w)).toArray();
      
      int capacity = 100;
      int shortageCost = 100;
      
      NormalDist[] weights = IntStream.iterate(0, i -> i + 1).limit(expectedWeights.length)
                                      .mapToObj(i -> new NormalDist(expectedWeights[i], standardDeviationWeights[i]))
                                      .toArray(NormalDist[]::new);
      
      SKPNormal instance = new SKPNormal(expectedValues, weights, capacity, shortageCost);
      
      int partitions = 100;
      SKPNormalMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new SKPNormalMILP(instance, partitions);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      
      SimulateNormal sim = new SimulateNormal(instance, seed);
      double milpSolutionValue = milp.getSolutionValue();
      double milpLinearizationError = milp.getMaxLinearizationError();
      int nbSamples = 1000;
      double simSolutionValue = sim.simulate(knapsack, nbSamples);
      
      System.out.println("MILP: "+milpSolutionValue);
      System.out.println("MILP max linearization error: "+milpLinearizationError);
      System.out.println("Simulation: "+simSolutionValue);
      System.out.println("Linearization gap (%): "+100*(simSolutionValue-milpSolutionValue)/simSolutionValue);
   }
}
