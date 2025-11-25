package skp.sim;

import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import ilog.concert.IloException;

import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.KP;
import skp.instance.SKPNormal;
import skp.milp.KPMILP;
import skp.milp.PWAPPROXIMATION;
import skp.milp.SKPNormalMILP;
import skp.sim.instance.SKPNormalRecedingSolvedInstance;
import skp.utilities.gson.GSONUtility;
import skp.utilities.probability.SampleFactory;

public class SimulateNormalReceding extends Simulate {
   SKPNormal instance;
   int partitions;
   double s = 1e-2; // sqrt approximation step
   
   public SimulateNormalReceding(SKPNormal instance, int partitions, double s) {
      super();
      this.instance = instance;
      this.partitions = partitions;
      this.s = s;
   }
   
   public SKPNormalRecedingSolvedInstance solve(int simulationRuns) {
      
      SimulateNormalReceding sim = new SimulateNormalReceding(instance, partitions, s);
      
      double[] realisations = sim.simulate(simulationRuns, partitions, s);
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
      
      SKPNormalRecedingSolvedInstance solvedInstance = new SKPNormalRecedingSolvedInstance(
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
   
   private int simulateOneItem(int t, double[] realizations, double remainingCapacity, int partitions, double s) {
      double[] previousRealizations = new double[t];
      System.arraycopy(realizations, 0, previousRealizations, 0, t);
      
      double[] expWeight = Arrays.stream(instance.getWeights()).mapToDouble(w -> w.getMean()).toArray();
      double[] stdWeight = Arrays.stream(instance.getWeights()).mapToDouble(w -> w.getSigma()).toArray();
      
      int Nbperiods = expWeight.length-t;
      
      double[] shortExpWeight = new double[Nbperiods];
      System.arraycopy(expWeight, t, shortExpWeight, 0, Nbperiods);
      
      double[] shortStdWeight = new double[Nbperiods];
      System.arraycopy(stdWeight, t, shortStdWeight, 0, Nbperiods);
      
      double[] shortExpValues = new double[Nbperiods];
      System.arraycopy(instance.getExpectedValues(), t, shortExpValues, 0, Nbperiods);
      
      SKPNormal reducedInstance = new SKPNormal(shortExpValues, shortExpWeight, shortStdWeight, remainingCapacity, instance.getShortageCost());
      
      SKPNormalMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new SKPNormalMILP(reducedInstance, partitions, s, PWAPPROXIMATION.EDMUNDSON_MADANSKI);
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
   
   private double simulateOneRun(double[] realizations, int partitions, double s) {
      double knapsackValue = 0;
      double remainingCapacity = instance.getCapacity();
      for(int i = 0; i < realizations.length; i++) {
         if(simulateOneItem(i, realizations, remainingCapacity, partitions, s) == 1) {
            remainingCapacity -= realizations[i];
            knapsackValue += this.instance.getExpectedValues()[i];
         }
      }
      knapsackValue -= Math.max(-remainingCapacity*instance.getShortageCost(), 0);
      return knapsackValue;
   }
   
   double[] simulate(int simulationRuns, int partitions, double s) {
      double[][] sampleMatrix = sampleWeights(simulationRuns);
      double[] knapsackValues = Arrays.stream(sampleMatrix)
                                      .parallel()
                                      .mapToDouble(r -> simulateOneRun(r, partitions, s))
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

      SKPNormal instance = SKPNormal.getTestInstance();
      
      int partitions = 10;
      double s = 1e-2;
      int simulationRuns = 100;
      
      SimulateNormalReceding sim = new SimulateNormalReceding(instance, partitions, s);
      System.out.println(GSONUtility.<SKPNormalRecedingSolvedInstance>printInstanceAsJSON(sim.solve(simulationRuns)));
   }
}
