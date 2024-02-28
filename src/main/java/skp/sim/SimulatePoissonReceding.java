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

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import ilog.concert.IloException;

import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.KP;
import skp.instance.SKPPoisson;
import skp.milp.KPMILP;
import skp.milp.SKPPoissonMILP;
import skp.sim.instance.SKPPoissonRecedingSolvedInstance;

import umontreal.ssj.randvar.UniformGen;


public class SimulatePoissonReceding extends Simulate {
   
   SKPPoisson instance;
   int partitions;
   
   public SimulatePoissonReceding(SKPPoisson instance, int partitions) {
      super();
      this.instance = instance;
      this.partitions = partitions;
   }
   
   public SKPPoissonRecedingSolvedInstance solve(int linearisationSamples, int simulationRuns) {
      
      SimulatePoissonReceding sim = new SimulatePoissonReceding(instance, partitions);
      
      double[] realisations = sim.simulate(simulationRuns, partitions, linearisationSamples, simulationRuns);
      Mean m = new Mean();
      double simSolutionMean = m.evaluate(realisations);
      StandardDeviation std = new StandardDeviation();
      double simSolutionStd = std.evaluate(realisations);
      double simEVwPI = sim.simulateEVwPI(simulationRuns);
      double EVP = sim.computeEVP();
      double EVwPI_obj_2_n = sim.simulateEVwPI_obj_2_n(simulationRuns);
      
      System.out.println("Simulation (mean): "+simSolutionMean);
      System.out.println("Simulation (std): "+simSolutionStd);
      System.out.println("EVwPI: "+simEVwPI);
      System.out.println("EVP: "+EVP);
      System.out.println("EVwPI on items 2-n: "+EVwPI_obj_2_n);
      
      SKPPoissonRecedingSolvedInstance solvedInstance = new SKPPoissonRecedingSolvedInstance(
            this.instance,
            simSolutionMean,
            simSolutionStd,
            simulationRuns,
            EVP,
            simEVwPI,
            EVwPI_obj_2_n,
            partitions,
            PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples()
            );
      
      return solvedInstance;
   }
   
   private int simulateOneItem(int t, double[] realizations, int remainingCapacity, int partitions, int linearisationSamples, int simulationRuns) {
      double[] previousRealizations = new double[t];
      System.arraycopy(realizations, 0, previousRealizations, 0, t);

      double[] expWeight = Arrays.stream(instance.getWeights()).mapToDouble(w -> w.getMean()).toArray();

      int Nbperiods = expWeight.length-t;

      double[] shortExpWeight = new double[Nbperiods];
      System.arraycopy(expWeight, t, shortExpWeight, 0, Nbperiods);
      
      double[] shortExpValues = new double[Nbperiods];
      System.arraycopy(instance.getExpectedValuesPerUnit(), t, shortExpValues, 0, Nbperiods);

      SKPPoisson reducedInstance = new SKPPoisson(shortExpValues, shortExpWeight, remainingCapacity, instance.getShortageCost());

      SKPPoissonMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new SKPPoissonMILP(reducedInstance, partitions, linearisationSamples);
         milp.solve(simulationRuns);
         knapsack = milp.getOptimalKnapsack();
         //System.out.println("Knapsack: "+Arrays.toString(knapsack));
         System.out.print(".");
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }

      return knapsack[0];
   }

   private double simulateOneRun(double[] realizations, int partitions, int linearisationSamples, int simulationRuns) {
      double knapsackValue = 0;
      int remainingCapacity = instance.getCapacity();
      for(int i = 0; i < realizations.length; i++) {
         if(simulateOneItem(i, realizations, remainingCapacity, partitions, linearisationSamples, simulationRuns) == 1) {
            remainingCapacity -= realizations[i];
            knapsackValue += this.instance.getExpectedValuesPerUnit()[i]*realizations[i];
         }
      }
      knapsackValue -= Math.max(-remainingCapacity*instance.getShortageCost(), 0);
      return knapsackValue;
   }
   
   double[] simulate(int nbSamples, int partitions, int linearisationSamples, int simulationRuns) {
      double[][] sampleMatrix = sampleWeights(nbSamples);
      double[] knapsackValues = Arrays.stream(sampleMatrix)
                                      .parallel()
                                      .mapToDouble(r -> simulateOneRun(r, partitions, linearisationSamples, simulationRuns))
                                      .peek(r -> System.out.println("Simulation run completed: "+r))
                                      .toArray();
      return knapsackValues;
   }
   
   private double simulateOneRunEVwPI(double[] realizations) {
      KP instanceEVwPI = new KP(instance.getExpectedValuesPerUnit(), realizations, instance.getCapacity(), instance.getShortageCost());
      KPMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new KPMILP(instanceEVwPI);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      return milp.getSolutionValue();
   }
   
   double simulateEVwPI(int nbSamples) {
      double[][] sampleMatrix = sampleWeights(nbSamples);
      double knapsackValue = Arrays.stream(sampleMatrix).parallel().mapToDouble(r -> simulateOneRunEVwPI(r)).sum()/nbSamples;
      return knapsackValue;
   }
   
   double computeEVP() {
      KP instanceEVP = new KP(instance.getExpectedValuesPerUnit(), Arrays.stream(instance.getWeights()).mapToDouble(w -> w.getMean()).toArray(), instance.getCapacity(), instance.getShortageCost());
      KPMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new KPMILP(instanceEVP);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      return milp.getSolutionValue();
   }
   
   private double simulateOneRunEVwPI(double[] realizations, int selectFirstObject) {
      KP instanceEVPI = new KP(instance.getExpectedValuesPerUnit(), realizations, instance.getCapacity(), instance.getShortageCost());
      KPMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new KPMILP(instanceEVPI, selectFirstObject);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      return milp.getSolutionValue();
   }
   
   double simulateEVwPI_obj_2_n(int nbSamples) {
      double[][] sampleMatrix = sampleWeights(nbSamples);
      double knapsackValueX1Eq0 = Arrays.stream(sampleMatrix).parallel().mapToDouble(r -> simulateOneRunEVwPI(r, 0)).sum()/nbSamples;
      double knapsackValueX1Eq1 = Arrays.stream(sampleMatrix).parallel().mapToDouble(r -> simulateOneRunEVwPI(r, 1)).sum()/nbSamples;
      return Math.max(knapsackValueX1Eq0, knapsackValueX1Eq1);
   }
   
   private double[][] sampleWeights(int nbSamples) {
      this.randGenerator.resetStartStream();
      double[][] sampleMatrix = new double[nbSamples][instance.getWeights().length];
      for(int i = 0; i < sampleMatrix.length; i++){
         for(int j = 0; j < sampleMatrix[i].length; j++){
            sampleMatrix[i][j] = instance.getWeights()[j].inverseF(UniformGen.nextDouble(this.randGenerator, 0, 1));
         }
      }
      return sampleMatrix;
   }
   
   public static void main(String args[]) {

      SKPPoisson instance = SKPPoisson.getTestInstance();
      
      int partitions = 10;
      int linearisationSamples = 1000;
      int simulationRuns = 1000;
      int nbSamples = 100;
      
      SimulatePoissonReceding sim = new SimulatePoissonReceding(instance, partitions);
      
      double[] realisations = sim.simulate(nbSamples, partitions, linearisationSamples, simulationRuns);
      Mean m = new Mean();
      double simSolutionMean = m.evaluate(realisations);
      StandardDeviation std = new StandardDeviation();
      double simSolutionStd = std.evaluate(realisations);
      double simEVwPI = sim.simulateEVwPI(nbSamples);
      double EVP = sim.computeEVP();
      double EVwPI_obj_2_n = sim.simulateEVwPI_obj_2_n(nbSamples);
      
      System.out.println("Simulation (mean): "+simSolutionMean);
      System.out.println("Simulation (std): "+simSolutionStd);
      System.out.println("Simulation (CI): ["+(simSolutionMean-1.96*simSolutionStd/Math.sqrt(nbSamples))+","+
                                              (simSolutionMean+1.96*simSolutionStd/Math.sqrt(nbSamples))+"]");
      System.out.println("EVwPI: "+simEVwPI);
      System.out.println("EVP: "+EVP);
      System.out.println("EVwPI on items 2-n: "+EVwPI_obj_2_n);
   }
}
