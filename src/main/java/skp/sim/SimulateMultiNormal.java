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
import skp.instance.SKPMultiNormal;
import skp.milp.SKPMultiNormalMILP;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalCholeskyGen;
import umontreal.ssj.randvarmulti.MultinormalGen;
import umontreal.ssj.randvarmulti.MultinormalPCAGen;
import umontreal.ssj.rng.MRG32k3aL;

public class SimulateMultiNormal {
   private static final long[] seed = {1,2,3,4,5,6};
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
         if(knapsack[i] == 1) knapsackValue += this.instance.getExpectedValuesPerUnit()[i]*this.instance.getWeights().getMean()[i]; 
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
      //MultinormalGen gen = new MultinormalCholeskyGen(standardNormal, mu, sigma);
      MultinormalGen gen = new MultinormalPCAGen(standardNormal, mu, sigma);
      double[][] points = new double[nbSamples][this.instance.getItems()];
      for(int i = 0; i < nbSamples; i++)
         gen.nextPoint(points[i]);
      return points;
   }
   
   
   
   public static void main(String args[]) {
      
      SKPMultiNormal instance = SKPMultiNormal.getTestInstance();
      
      int partitions = 20;
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
      double milpLinearizationError = milp.getMaxLinearizationError();
      int nbSamples = 10000;
      double simSolutionValue = sim.simulate(knapsack, nbSamples);
      
      System.out.println("MILP: "+milpSolutionValue);
      System.out.println("MILP max linearization error: "+100*milpLinearizationError/simSolutionValue);
      System.out.println("Simulation: "+simSolutionValue);
      System.out.println("Linearization gap (%): "+100*(simSolutionValue-milpSolutionValue)/simSolutionValue);
   }

}
