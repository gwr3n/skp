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
import java.util.Arrays;

import ilog.concert.IloException;
import skp.instance.SKPGenericDistribution;
import skp.milp.SKPGenericDistributionCuts;
import skp.milp.instance.SKPGenericDistributionBandBSolvedInstance;
import skp.milp.instance.SKPGenericDistributionCutsSolvedInstance;
import skp.sdp.DSKPGenericDistribution;
import skp.sdp.instance.DSKPGenericDistributionSolvedInstance;
import skp.utilities.gson.GSONUtility;

import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.randvar.RandomVariateGen;

public class SKPGenericDistributionBatch extends SKPBatch {
   
   public static void main(String args[]) {
      
      int[] instanceSize = {25, 50, 100, 500};
      double[] coeff_of_var  = {0.1, 0.2};
      INSTANCE_TYPE[] instanceType = {
            INSTANCE_TYPE.P05_UNCORRELATED,
            INSTANCE_TYPE.P05_WEAKLY_CORRELATED,
            INSTANCE_TYPE.P05_STRONGLY_CORRELATED,
            INSTANCE_TYPE.P05_INVERSE_STRONGLY_CORRELATED,
            INSTANCE_TYPE.P05_ALMOST_STRONGLY_CORRELATED,
            INSTANCE_TYPE.P05_SUBSET_SUM,
            INSTANCE_TYPE.P05_UNCORRELATED_SIMILAR_WEIGHTS,
            INSTANCE_TYPE.P05_PROFIT_CEILING,
            INSTANCE_TYPE.P05_CIRCLE_INSTANCES};
      
      for(INSTANCE_TYPE t: instanceType) {
         for(int size : instanceSize) {
            for(double cv : coeff_of_var) {
               File folder = new File("batch/"+t.toString()+"/"+size+"/"+cv);
               if (!folder.exists()) {
                  folder.mkdirs();
               }
               
               String batchFileName = "batch/"+t.toString()+"/"+size+"/"+cv+"/generic_distribution_instances.json";
               
               SKPGenericDistribution[] instances = generateInstances(batchFileName, t, size, cv);
               
               int linearizationSamples = 500;
               int simulationRuns = 500;   
               int maxCuts = 1000;
               try {
                  solveMILP(instances, linearizationSamples, maxCuts, simulationRuns, "batch/"+t.toString()+"/"+size+"/"+cv);
               } catch (IloException e) {
                  e.printStackTrace();
               }
               if(size == instanceSize[0])
                  solveDSKP(instances, "batch/"+t.toString()+"/"+size+"/"+cv);
            }
         }
      }
      
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
   
   private static void bubbleSort(double[] expectedValues, double[] expectedWeights) {
      if (expectedValues.length != expectedWeights.length) {
         throw new IllegalArgumentException("Arrays must be of the same length");
      }

      boolean swapped;
      do {
         swapped = false;
         for (int i = 0; i < expectedValues.length - 1; i++) {
            double ratioA = expectedValues[i] / expectedWeights[i];
            double ratioB = expectedValues[i + 1] / expectedWeights[i + 1];

            if (ratioA < ratioB) {
               // Swap values
               double tempValue = expectedValues[i];
               expectedValues[i] = expectedValues[i + 1];
               expectedValues[i + 1] = tempValue;

               // Swap weights
               double tempWeight = expectedWeights[i];
               expectedWeights[i] = expectedWeights[i + 1];
               expectedWeights[i + 1] = tempWeight;

               swapped = true;
            }
         }
      } while (swapped);
   }
   
   /**
    * Generate a batch of instances
    */

   protected static SKPGenericDistribution[] generateInstances(String batchFileName, INSTANCE_TYPE type, int instanceSize, double cv) {
      int H = 10;
      int R = 100;
      double shortageCost = 10;
      
      switch(type) {
         /*case GAMMA: {
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
         }*/
         case P05_UNCORRELATED: {
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               
               bubbleSort(expectedValues, expectedWeights);
               
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
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double U = Math.max(1.0, Arrays.stream(expectedWeights).map(v -> v - R/10.0).max().getAsDouble());
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(U,U+2*R/10.0)).nextArrayOfDouble(instanceSize);
               
               bubbleSort(expectedValues, expectedWeights);
               
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
         case P05_STRONGLY_CORRELATED: {
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> v + R/10.0).toArray();
               
               bubbleSort(expectedValues, expectedWeights);
               
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
         case P05_INVERSE_STRONGLY_CORRELATED: {
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedWeights = Arrays.stream(expectedValues).map(v -> v + R/10.0).toArray();
               
               bubbleSort(expectedValues, expectedWeights);
               
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
         case P05_ALMOST_STRONGLY_CORRELATED: {
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> new RandomVariateGen(randGenerator, new UniformDist(v + R/10.0 - R/500.0, v + R/10.0 + R/500.0)).nextDouble()).toArray();
               
               bubbleSort(expectedValues, expectedWeights);
               
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
         case P05_SUBSET_SUM: {
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.copyOf(expectedWeights, instanceSize);
               
               bubbleSort(expectedValues, expectedWeights);
               
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
         case P05_UNCORRELATED_SIMILAR_WEIGHTS: {
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(R,R+10)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               
               bubbleSort(expectedValues, expectedWeights);
               
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
         case P05_PROFIT_CEILING: {
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double d = 3;
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> d*Math.ceil(v/d)).toArray();
               
               bubbleSort(expectedValues, expectedWeights);
               
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
         case P05_CIRCLE_INSTANCES: {
            
            SKPGenericDistribution[] batch = new SKPGenericDistribution[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> 2*Math.sqrt(4*R*R - Math.pow(v - 2*R,2))/3).toArray();
               
               bubbleSort(expectedValues, expectedWeights);
               
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
   
   public static void solveMILP(SKPGenericDistribution[] batch, int linearizationSamples, int maxCuts, int simulationRuns, String folder) throws IloException {
      // SKPGenericDistribution[] batch = retrieveBatch(fileName); // Batch cannot be retrieved because Distribution[] is not Serializable
      
      String fileNameSolved = folder+"/solved_generic_distribution_instances_MILP.json";
      SKPGenericDistributionCutsSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, linearizationSamples, maxCuts, simulationRuns);
      
      // solvedBatch = retrieveSolvedBatchMILP(fileNameSolved); // Batch cannot be retrieved because Distribution[] is not Serializable
      System.out.println(GSONUtility.<SKPGenericDistributionCutsSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = folder+"/solved_generic_distribution_instances_MILP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPGenericDistributionCutsSolvedInstance[] solveBatchMILP(SKPGenericDistribution[] instances, String fileName, int linearizationSamples, int maxCuts, int simulationRuns) throws IloException {
      /*
       * Sequential
       * 
      ArrayList<SKPGenericDistributionBandBSolvedInstance>solved = new ArrayList<SKPGenericDistributionBandBSolvedInstance>();
      for(SKPGenericDistribution instance : instances) {
         solved.add(new SKPGenericDistributionBandB(instance, linearizationSamples, simulationRuns).solve());
         
      }
      GSONUtility.<SKPGenericDistributionBandBSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPGenericDistributionBandBSolvedInstance[solved.size()]), fileName);
      return solved.toArray(new SKPGenericDistributionBandBSolvedInstance[solved.size()]);
      */
      
      /*
       * Parallel
       */
      SKPGenericDistributionCutsSolvedInstance[] solved = Arrays.stream(instances)
                                                                .parallel()
                                                                .map(instance -> {
         try {
            return new SKPGenericDistributionCuts(instance, linearizationSamples, maxCuts, simulationRuns).solve();
         } catch (IloException e) {
            e.printStackTrace();
            return null;
         }
      }).toArray(SKPGenericDistributionCutsSolvedInstance[]::new);
      GSONUtility.<SKPGenericDistributionCutsSolvedInstance[]>saveInstanceToJSON(solved, fileName);
      return solved;
   }
   
   /*
    * Batch cannot be retrieved because Distribution[] is not Serializable
    *
   private static SKPGenericDistributionMILPSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPGenericDistributionMILPSolvedInstance[] solvedInstances = GSONUtility.<SKPGenericDistributionMILPSolvedInstance[]>retrieveJSONInstance(fileName, SKPGenericDistributionMILPSolvedInstance[].class);
      return solvedInstances;
   }*/
   
   @SuppressWarnings("unused")
   private static void storeSolvedBatchToCSV(SKPGenericDistributionCutsSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, milpOptimalityGap, cuts, "
            + "piecewiseSamples, milpMaxLinearizationError, simulatedLinearizationError,"
            + "cplexSolutionTimeMs, simplexIterations, exploredNodes\n";
      String body = "";
      
      for(SKPGenericDistributionCutsSolvedInstance s : instances) {
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
   
   @SuppressWarnings("unused")
   private static void storeSolvedBatchToCSV(SKPGenericDistributionBandBSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, solutionTimeMs, exploredNodes, optGap\n";
      String body = "";
      
      for(SKPGenericDistributionBandBSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).map(d -> d.toString()).toArray()).replace(",", "\t")+ ", " +
                 s.instance.getCapacity()+ ", " +
                 s.instance.getShortageCost()+ ", " +
                 Arrays.toString(s.optimalKnapsack).replace(",", "\t")+ ", " +
                 s.simulatedSolutionValue + ", " +
                 s.simulationRuns + ", " +
                 s.solutionTimeMs + ", " +
                 s.exploredNodes + ", " +
                 s.optGap +"\n";
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
   
   public static void solveDSKP(SKPGenericDistribution[] batch, String folder) {
      // SKPNormal[] batch = retrieveBatch(fileName); // Batch cannot be retrieved because Distribution[] is not Serializable
      
      String fileNameSolved = folder+"/solved_generic_distribution_instances_DSKP.json";
      DSKPGenericDistributionSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      // solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved); // Batch cannot be retrieved because Distribution[] is not Serializable
      System.out.println(GSONUtility.<DSKPGenericDistributionSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = folder+"/solved_generic_distribution_instances_DSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static DSKPGenericDistributionSolvedInstance[] solveBatchDSKP(SKPGenericDistribution[] instances, String fileName) {
      /*
       * Sequential
       * 
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPGenericDistributionSolvedInstance>solved = new ArrayList<DSKPGenericDistributionSolvedInstance>();
      for(SKPGenericDistribution instance : instances) {
         solved.add(new DSKPGenericDistribution(instance, truncationQuantile).solve());
         System.out.println("Solved DSKP instance number "+solved.size());
         GSONUtility.<DSKPGenericDistributionSolvedInstance[]>saveInstanceToJSON(solved.toArray(new DSKPGenericDistributionSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPGenericDistributionSolvedInstance[solved.size()]);*/
      
      /*
       * Parallel
       */
      double truncationQuantile = 0.999999999999999;
      DSKPGenericDistributionSolvedInstance[] solved = Arrays.stream(instances)
                                                             .parallel()
                                                             .map(instance -> new DSKPGenericDistribution(instance, truncationQuantile).solve())
                                                             .toArray(DSKPGenericDistributionSolvedInstance[]::new);
      GSONUtility.<DSKPGenericDistributionSolvedInstance[]>saveInstanceToJSON(solved, fileName);
      return solved;
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
