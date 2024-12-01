package skp.batch;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import ilog.concert.IloException;

import skp.instance.SKPGenericDistribution;
import skp.sim.SimulateGenericDistributionReceding;
import skp.sim.instance.SKPGenericDistributionRecedingSolvedInstance;
import skp.utilities.gson.GSONUtility;

public class SKPGenericDistributionRecedingBatch extends SKPGenericDistributionBatch{
   
   public static void main(String args[]) {
      
      int[] instanceSize = {25};
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
               
               int linearizationSamples = 1000;
               int simulationRunsRH = 100; 
               int simulationRuns = 10000;
               int maxCuts = 1000;
               try {
                  solveMILP(instances, simulationRunsRH, linearizationSamples, maxCuts, simulationRuns, "batch/"+t.toString()+"/"+size+"/"+cv);
               } catch (IloException e) {
                  e.printStackTrace();
               }
            }
         }
      }
      
      
   }

   /*
    * MILP receding
    */ 
   
   public static void solveMILP(SKPGenericDistribution[] batch, int simulationRunsRH, int linearizationSamples, int maxCuts, int simulationRuns, String folder) throws IloException {
      // SKPGenericDistribution[] batch = retrieveBatch(fileName); // Batch cannot be retrieved because Distribution[] is not Serializable
      
      String fileNameSolved = folder+"/solved_generic_distribution_instances_MILP_receding.json";
      SKPGenericDistributionRecedingSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, linearizationSamples, simulationRunsRH, maxCuts, simulationRuns);
      
      //solvedBatch = retrieveSolvedBatchMILP(fileNameSolved); // Batch cannot be retrieved because Distribution[] is not Serializable
      System.out.println(GSONUtility.<SKPGenericDistributionRecedingSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = folder+"/solved_generic_distribution_instances_MILP_receding.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPGenericDistributionRecedingSolvedInstance[] solveBatchMILP(SKPGenericDistribution[] instances, String fileName, int linearizationSamples, int simulationRunsRH, int maxCuts, int simulationRuns) throws IloException {
      ArrayList<SKPGenericDistributionRecedingSolvedInstance>solved = new ArrayList<SKPGenericDistributionRecedingSolvedInstance>();
      for(SKPGenericDistribution instance : instances) {
         solved.add(new SimulateGenericDistributionReceding(instance).solve(simulationRunsRH, linearizationSamples, maxCuts, simulationRuns));
         System.out.println("Solved instance number "+solved.size());
         GSONUtility.<SKPGenericDistributionRecedingSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPGenericDistributionRecedingSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPGenericDistributionRecedingSolvedInstance[solved.size()]);
   }
   
   private static void storeSolvedBatchToCSV(SKPGenericDistributionRecedingSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights, "
            + "capacity, shortageCost, simulatedSolutionMean, simulatedSolutionStd, "
            + "simulationRuns, "
            + "piecewiseSamples, evp, evwpi,"
            + "evwpi_obj_2_n\n";
      String body = "";
      
      for(SKPGenericDistributionRecedingSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).map(d -> d.toString()).toArray()).replace(",", "\t")+ ", " +
                 s.instance.getCapacity()+ ", " +
                 s.instance.getShortageCost()+ ", " +
                 s.simulatedSolutionMean + ", " +
                 s.simulatedSolutionStd + ", " +
                 s.simulationRuns + ", " +
                 s.linearisationSamples + ", " +
                 s.evp + ", " +
                 s.evwpi + ", " +
                 s.evwpi_obj_2_n +"\n";
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
}
