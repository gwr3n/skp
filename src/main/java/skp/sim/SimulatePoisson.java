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
import skp.SKPPoisson;
import skp.milp.SKPPoissonMILP;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.rng.MRG32k3aL;

public class SimulatePoisson {
   SKPPoisson instance;
   private MRG32k3aL randGenerator;
   
   public SimulatePoisson(SKPPoisson instance, long[] seed) {
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
      long[] seed = {1,2,3,4,5,6};
      
      double[] expectedValues = {100,100,100,100,100};
      double[] expectedWeights = {10,10,40,30,20};
      
      int capacity = 100;
      int shortageCost = 100;
      
      PoissonDist[] weights = IntStream.iterate(0, i -> i + 1).limit(expectedWeights.length)
                                       .mapToObj(i -> new PoissonDist(expectedWeights[i]))
                                       .toArray(PoissonDist[]::new);
      
      SKPPoisson instance = new SKPPoisson(expectedValues, weights, capacity, shortageCost);
      
      int partitions = 10;
      int linearizationSamples = 10000;
      SKPPoissonMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new SKPPoissonMILP(instance, partitions, linearizationSamples);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      
      int simulationSamples = 100000;
      SimulatePoisson sim = new SimulatePoisson(instance, seed);
      double milpSolutionValue = milp.getSolutionValue();
      double milpLinearizationError = milp.getMaxLinearizationError();
      double simSolutionValue = sim.simulate(knapsack, simulationSamples);
      
      System.out.println("MILP: "+milpSolutionValue);
      System.out.println("MILP max linearization error: "+milpLinearizationError);
      System.out.println("Simulation: "+simSolutionValue);
      System.out.println("Linearization gap (%): "+100*(simSolutionValue-milpSolutionValue)/simSolutionValue);
   }
}
