/**
 * To run from Mac OS
 * 
 * -Djava.library.path=/Applications/CPLEX_Studio2211/opl/bin/x86-64_osx/
 * 
 * Environment variable
 * 
 * DYLD_LIBRARY_PATH /Applications/CPLEX_Studio2211/opl/bin/x86-64_osx/
 * 
 * @author Roberto Rossi
 *
 */

package skp.batch;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

import ilog.concert.IloException;

import skp.instance.SKPGenericDistribution;
import skp.milp.SKPGenericDistributionMILP;
import skp.milp.instance.SKPGenericDistributionMILPSolvedInstance;
import skp.sdp.DSKPGenericDistribution;
import skp.sdp.instance.DSKPGenericDistributionSolvedInstance;
import skp.utilities.gson.GSONUtility;

import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.probdist.UniformIntDist;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.randvar.RandomVariateGenInt;

public class SKPGenericDistributionBatch extends SKPBatch {
   
   public static void main(String args[]) {
      File folder = new File("batch");
      if (!folder.exists()) {
        folder.mkdir();
      } 
      
      String batchFileName = "batch/generic_distribution_instances.json";
      
      SKPGenericDistribution[] instances = generateInstances(batchFileName, INSTANCE_TYPE.P05_UNCORRELATED);
      
      int linearizationSamples = 100000;
      int simulationRuns = 100000;   
      int maxCuts = 1000;
      try {
         solveMILP(instances, linearizationSamples, simulationRuns, maxCuts);
      } catch (IloException e) {
         e.printStackTrace();
      }
      solveDSKP(instances);
   }
   
   enum INSTANCE_TYPE {
      GAMMA,
      P05_UNCORRELATED,
      P05_WEAKLY_CORRELATED,
      P05_STRONGLY_CORRELATED,
      P05_INVERSE_STRONGLY_CORRELATED,
      P05_ALMOST_STRONGLY_CORRELATED,
      P05_SUBSET_SUM,
      P05_UNCORRELATED_SIMILAR_WEIGHTS,
      P05_PROFIT_CEILING,
      P05_CIRCLE_INSTANCES
   };
   
   /**
    * Generate a batch of instances
    */

   protected static SKPGenericDistribution[] generateInstances(String batchFileName, INSTANCE_TYPE type) {
      
      switch(type) {
         case GAMMA: {
            int instances = 10;
            int instanceSize = 10;
            
            Distribution expectedValue = new UniformDist(2.75,275);
            Distribution expectedWeight = new UniformDist(15,70);
            double coefficientOfVariation = 0.2;
            DiscreteDistributionInt capacity = new UniformIntDist(100,200);
            Distribution shortageCost = new UniformDist(2,10);
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[instances];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            batch = IntStream.iterate(0, i -> i + 1)
                             .limit(instances)
                             .mapToObj(i -> new SKPGenericDistribution(
                                                                (new RandomVariateGen(randGenerator, expectedValue)).nextArrayOfDouble(instanceSize),
                                                                Arrays.stream((new RandomVariateGen(randGenerator, expectedWeight)).nextArrayOfDouble(instanceSize)).mapToObj(w -> new GammaDist(1/Math.pow(coefficientOfVariation,2), 1/(w*Math.pow(coefficientOfVariation,2)))).toArray(GammaDist[]::new),
                                                                (new RandomVariateGenInt(randGenerator, capacity)).nextInt(),
                                                                (new RandomVariateGen(randGenerator, shortageCost)).nextDouble()))
                             .toArray(SKPGenericDistribution[]::new);
            
            return batch;
         }
         case P05_UNCORRELATED: {
            int H = 10;
            int instanceSize = 10;
            int R = 100;
            double cv = 0.2;
            double shortageCost = 10;
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPGenericDistribution(
                     expectedValues,
                     Arrays.stream(expectedWeights).mapToObj(w -> new GammaDist(1/Math.pow(cv,2), 1/(w*Math.pow(cv,2)))).toArray(GammaDist[]::new),
                     capacity,
                     shortageCost
                     );
            }
            return batch;
         }
         case P05_WEAKLY_CORRELATED: {
            
         }
         case P05_STRONGLY_CORRELATED: {
            
         }
         case P05_INVERSE_STRONGLY_CORRELATED: {
            
         }
         case P05_ALMOST_STRONGLY_CORRELATED: {
            
         }
         case P05_SUBSET_SUM: {
            
         }
         case P05_UNCORRELATED_SIMILAR_WEIGHTS: {
            
         }
         case P05_PROFIT_CEILING: {
            
         }
         case P05_CIRCLE_INSTANCES: {
            
         }
         default:
            return null;
      }
   }
   
   /*
    * Batch cannot be retrieved because Distribution[] is not Serializable
    *
   public static SKPGenericDistribution[] retrieveBatch(String fileName) {
      SKPGenericDistribution[] instances = GSONUtility.<SKPGenericDistribution[]>retrieveJSONInstance(fileName, SKPGenericDistribution[].class);
      return instances;
   }*/
   
   public static void solveMILP(SKPGenericDistribution[] batch, int linearizationSamples, int simulationRuns, int maxCuts) throws IloException {
      // SKPGenericDistribution[] batch = retrieveBatch(fileName); // Batch cannot be retrieved because Distribution[] is not Serializable
      
      String fileNameSolved = "batch/solved_generic_distribution_instances_MILP.json";
      SKPGenericDistributionMILPSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, linearizationSamples, simulationRuns, maxCuts);
      
      // solvedBatch = retrieveSolvedBatchMILP(fileNameSolved); // Batch cannot be retrieved because Distribution[] is not Serializable
      System.out.println(GSONUtility.<SKPGenericDistributionMILPSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "batch/solved_generic_distribution_instances_MILP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPGenericDistributionMILPSolvedInstance[] solveBatchMILP(SKPGenericDistribution[] instances, String fileName, int linearizationSamples, int simulationRuns, int maxCuts) throws IloException {
      ArrayList<SKPGenericDistributionMILPSolvedInstance>solved = new ArrayList<SKPGenericDistributionMILPSolvedInstance>();
      for(SKPGenericDistribution instance : instances) {
         solved.add(new SKPGenericDistributionMILP(instance, linearizationSamples, maxCuts).solve(simulationRuns));
         GSONUtility.<SKPGenericDistributionMILPSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPGenericDistributionMILPSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPGenericDistributionMILPSolvedInstance[solved.size()]);
   }
   
   /*
    * Batch cannot be retrieved because Distribution[] is not Serializable
    *
   private static SKPGenericDistributionMILPSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPGenericDistributionMILPSolvedInstance[] solvedInstances = GSONUtility.<SKPGenericDistributionMILPSolvedInstance[]>retrieveJSONInstance(fileName, SKPGenericDistributionMILPSolvedInstance[].class);
      return solvedInstances;
   }*/
   
   private static void storeSolvedBatchToCSV(SKPGenericDistributionMILPSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, milpOptimalityGap, cuts, "
            + "piecewiseSamples, milpMaxLinearizationError, simulatedLinearizationError,"
            + "cplexSolutionTimeMs, simplexIterations, exploredNodes\n";
      String body = "";
      
      for(SKPGenericDistributionMILPSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).map(d -> d.toString()).toArray()).replace(",", "\t")+ ", " +
                 s.instance.getCapacity()+ ", " +
                 s.instance.getShortageCost()+ ", " +
                 Arrays.toString(s.optimalKnapsack).replace(",", "\t")+ ", " +
                 s.simulatedSolutionValue + ", " +
                 s.simulationRuns + ", " +
                 s.milpSolutionValue + ", " +
                 s.milpOptimalityGap + ", " +
                 s.cuts + ", " +
                 s.piecewiseSamples  + ", " +
                 s.milpMaxLinearizationError + ", " +
                 s.simulatedLinearizationError + ", " +
                 s.cplexSolutionTimeMs + ", " +
                 s.simplexIterations + ", " +
                 s.exploredNodes +"\n";
      }
      PrintWriter pw;
      try {
         pw = new PrintWriter(new File(fileName));
         pw.print(header+body);
         pw.close();
      } catch (FileNotFoundException e) {
         e.printStackTrace();
      }
   }
   
   /*
    * DSKP
    */
   
   public static void solveDSKP(SKPGenericDistribution[] batch) {
      // SKPNormal[] batch = retrieveBatch(fileName); // Batch cannot be retrieved because Distribution[] is not Serializable
      
      String fileNameSolved = "batch/solved_generic_distribution_instances_DSKP.json";
      DSKPGenericDistributionSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      // solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved); // Batch cannot be retrieved because Distribution[] is not Serializable
      System.out.println(GSONUtility.<DSKPGenericDistributionSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "batch/solved_generic_distribution_instances_DSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static DSKPGenericDistributionSolvedInstance[] solveBatchDSKP(SKPGenericDistribution[] instances, String fileName) {
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPGenericDistributionSolvedInstance>solved = new ArrayList<DSKPGenericDistributionSolvedInstance>();
      for(SKPGenericDistribution instance : instances) {
         solved.add(new DSKPGenericDistribution(instance, truncationQuantile).solve());
         GSONUtility.<DSKPGenericDistributionSolvedInstance[]>saveInstanceToJSON(solved.toArray(new DSKPGenericDistributionSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPGenericDistributionSolvedInstance[solved.size()]);
   }

   private static void storeSolvedBatchToCSV(DSKPGenericDistributionSolvedInstance[] instances, String fileName) {
      String header = "instanceID, expectedValues, expectedWeights, capacity, shortageCost, solutionValue, solutionTimeMs, statesExplored\n";
      String body = "";
      
      for(DSKPGenericDistributionSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).map(d -> d.toString()).toArray()).replace(",", "\t")+ ", " +
                 s.instance.getCapacity()+ ", " +
                 s.instance.getShortageCost()+ ", " +
                 s.solutionValue + ", " +
                 s.solutionTimeMs + ", " +
                 s.statesExplored +"\n";
      }
      PrintWriter pw;
      try {
         pw = new PrintWriter(new File(fileName));
         pw.print(header+body);
         pw.close();
      } catch (FileNotFoundException e) {
         e.printStackTrace();
      }
   }
   
   /*
    * Batch cannot be retrieved because Distribution[] is not Serializable
    *
   private static DSKPGenericDistributionSolvedInstance[] retrieveSolvedBatchDSKP(String fileName) {
      DSKPGenericDistributionSolvedInstance[] solvedInstances = GSONUtility.<DSKPGenericDistributionSolvedInstance[]>retrieveJSONInstance(fileName, DSKPGenericDistributionSolvedInstance[].class);
      return solvedInstances;
   }*/
   
}
