package skp.sim;

import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import ilog.concert.IloException;

import skp.instance.KP;
import skp.instance.SKPGenericDistribution;
import skp.milp.KPMILP;
import skp.milp.SKPGenericDistributionMILP;
import skp.sim.instance.SKPGenericDistributionRecedingSolvedInstance;
import skp.utilities.gson.GSONUtility;
import skp.utilities.probability.SampleFactory;

import umontreal.ssj.probdist.Distribution;

public class SimulateGenericDistributionReceding extends Simulate {
   SKPGenericDistribution instance;
   
   public SimulateGenericDistributionReceding(SKPGenericDistribution instance) {
      super();
      this.instance = instance;
   }
   
   public SKPGenericDistributionRecedingSolvedInstance solve(int simulationRuns, int linearisationSamples) {
      
      SimulateGenericDistributionReceding sim = new SimulateGenericDistributionReceding(instance);
      
      double[] realisations = sim.simulate(simulationRuns, linearisationSamples);
      Mean m = new Mean();
      double simSolutionMean = m.evaluate(realisations);
      StandardDeviation std = new StandardDeviation();
      double simSolutionStd = std.evaluate(realisations);
      double simEVwPI = sim.simulateEVwPI(simulationRuns);
      double EVP = sim.computeEVP();
      double EVwPI_obj_2_n = sim.simulateEVwPI_obj_2_n(simulationRuns);
      
      System.out.println("Simulation (mean): "+simSolutionMean);
      System.out.println("Simulation (std): "+simSolutionStd);
      System.out.println("Simulation (CI): ["+(simSolutionMean-1.96*simSolutionStd/Math.sqrt(simulationRuns))+","+
                                              (simSolutionMean+1.96*simSolutionStd/Math.sqrt(simulationRuns))+"]");
      System.out.println("EVwPI: "+simEVwPI);
      System.out.println("EVP: "+EVP);
      System.out.println("EVwPI on items 2-n: "+EVwPI_obj_2_n);
      
      SKPGenericDistributionRecedingSolvedInstance solvedInstance = new SKPGenericDistributionRecedingSolvedInstance(
            this.instance,
            simSolutionMean,
            simSolutionStd,
            simulationRuns,
            EVP,
            simEVwPI,
            EVwPI_obj_2_n,
            linearisationSamples
            );
      
      return solvedInstance;
   }
   
   private int simulateOneItem(int t, double[] realizations, double remainingCapacity, int linearisationSamples) {
      double[] previousRealizations = new double[t];
      System.arraycopy(realizations, 0, previousRealizations, 0, t);

      int Nbperiods = instance.getWeights().length-t;

      Distribution[] shortWeight = new Distribution[Nbperiods];
      System.arraycopy(instance.getWeights(), t, shortWeight, 0, Nbperiods);
      
      double[] shortExpValues = new double[Nbperiods];
      System.arraycopy(instance.getExpectedValues(), t, shortExpValues, 0, Nbperiods);

      SKPGenericDistribution reducedInstance = new SKPGenericDistribution(shortExpValues, shortWeight, remainingCapacity, instance.getShortageCost());

      SKPGenericDistributionMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new SKPGenericDistributionMILP(reducedInstance, linearisationSamples);
         milp.solve();
         knapsack = milp.getOptimalKnapsack();
         //System.out.println("Knapsack: "+Arrays.toString(knapsack));
         System.out.print(".");
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }

      return knapsack[0];
   }
   
   private double simulateOneRun(double[] realizations, int linearisationSamples) {
      double knapsackValue = 0;
      double remainingCapacity = instance.getCapacity();
      for(int i = 0; i < realizations.length; i++) {
         if(simulateOneItem(i, realizations, remainingCapacity, linearisationSamples) == 1) {
            remainingCapacity -= realizations[i];
            knapsackValue += this.instance.getExpectedValues()[i];
         }
      }
      knapsackValue -= Math.max(-remainingCapacity*instance.getShortageCost(), 0);
      return knapsackValue;
   }
   
   double[] simulate(int simulationRuns, int linearisationSamples) {
      double[][] sampleMatrix = sampleWeights(simulationRuns);
      double[] knapsackValues = Arrays.stream(sampleMatrix)
                                      .parallel()
                                      .mapToDouble(r -> simulateOneRun(r, linearisationSamples))
                                      .peek(r -> System.out.println("Simulation run completed: "+r))
                                      .toArray();
      return knapsackValues;
   }
   
   private double simulateOneRunEVwPI(double[] realizations) {
      KP instanceEVwPI = new KP(instance.getExpectedValues(), realizations, instance.getCapacity(), instance.getShortageCost());
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
      KP instanceEVP = new KP(instance.getExpectedValues(), Arrays.stream(instance.getWeights()).mapToDouble(w -> w.getMean()).toArray(), instance.getCapacity(), instance.getShortageCost());
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
      KP instanceEVPI = new KP(instance.getExpectedValues(), realizations, instance.getCapacity(), instance.getShortageCost());
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
      double[][] sampleMatrix;
      switch(Simulate.samplingStrategy) {
      case LHS:
         sampleMatrix = SampleFactory.getNextLHSample(instance.getWeights(), nbSamples, randGenerator);
         break;
      case SRS:
      default:
         sampleMatrix = SampleFactory.getNextSimpleRandomSample(instance.getWeights(), nbSamples, randGenerator);
      }
      return sampleMatrix;
   }
   
   public static void main(String args[]) {

      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstance();
      
      int linearisationSamples = 1000;
      int simulationRuns = 100;
      
      SimulateGenericDistributionReceding sim = new SimulateGenericDistributionReceding(instance);
      System.out.println(GSONUtility.<SKPGenericDistributionRecedingSolvedInstance>printInstanceAsJSON(sim.solve(simulationRuns, linearisationSamples)));
   }
}
