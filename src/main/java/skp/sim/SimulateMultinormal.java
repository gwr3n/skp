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
import skp.milp.SKPMultinormalMILP;
import skp.milp.instance.SKPMultinormalMILPSolvedInstance;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalGen;
import umontreal.ssj.randvarmulti.MultinormalPCAGen;

public class SimulateMultinormal extends Simulate {
   
   SKPMultinormal instance;
   
   public SimulateMultinormal(SKPMultinormal instance) {
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
   
   private double[][] sampleWeights(int nbSamples) {
      double[] mu = this.instance.getWeights().getMean();
      double[][] sigma = this.instance.getWeights().getCovariance();
      
      this.randGenerator.resetStartStream();
      NormalGen standardNormal = new NormalGen(this.randGenerator, 0, 1);
      MultinormalGen gen = new MultinormalPCAGen(standardNormal, mu, sigma);
      double[][] points = new double[nbSamples][this.instance.getItems()];
      for(int i = 0; i < nbSamples; i++)
         gen.nextPoint(points[i]);
      return points;
   }
   
   public static void main(String args[]) {
      
      SKPMultinormal instance = SKPMultinormal.getTestInstance();
      
      int partitions = 20;
      int simulationRuns = 10000;
      
      try {
         SKPMultinormalMILP milp = new SKPMultinormalMILP(instance, partitions);
         SKPMultinormalMILPSolvedInstance solved = milp.solve(simulationRuns);
         System.out.println(GSONUtility.<SKPMultinormalMILPSolvedInstance>printInstanceAsJSON(solved));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
   }

}
