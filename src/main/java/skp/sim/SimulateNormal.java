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
import skp.milp.SKPNormalMILP;
import skp.milp.instance.SKPNormalMILPSolvedInstance;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.UniformGen;

public class SimulateNormal extends Simulate {
   
   SKPNormal instance;
   
   public SimulateNormal(SKPNormal instance) {
      super();
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
      
      SKPNormal instance = SKPNormal.getTestInstance();
      
      int partitions = 10;
      int simulationSamples = 100000;
      
      try {
         SKPNormalMILP milp = new SKPNormalMILP(instance, partitions);
         SKPNormalMILPSolvedInstance solved = milp.solve(simulationSamples);
         System.out.println(GSONUtility.<SKPNormalMILPSolvedInstance>printInstanceAsJSON(solved));
      } catch (IloException e) {
         e.printStackTrace();
      }
   }
}
