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

import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;

public class SKPGenericDistributionRecedingBatch extends SKPGenericDistributionBatch{
   
   public static void main(String args[]) {
      File folder = new File("batch");
      if (!folder.exists()) {
        folder.mkdir();
      } 
      
      String batchFileName = "batch/generic_distribution_instances.json";
      
      SKPGenericDistribution[] instances = generateInstances(batchFileName);
      
      int linearizationSamples = 10000;
      int simulationRuns = 100; 
      try {
         solveMILP(instances, simulationRuns, linearizationSamples);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
     
   public SKPGenericDistributionRecedingBatch(
         Distribution expectedValuePerUnit,
         Distribution expectedWeight,
         Distribution coefficientOfVariation,
         DiscreteDistributionInt capacity,
         Distribution shortageCost) {
      super(expectedValuePerUnit, expectedWeight, coefficientOfVariation, capacity, shortageCost);
   }

   /*
    * MILP receding
    */ 
   
   public static void solveMILP(SKPGenericDistribution[] batch, int simulationRuns, int linearizationSamples) throws IloException {
      // SKPGenericDistribution[] batch = retrieveBatch(fileName); // Batch cannot be retrieved because Distribution[] is not Serializable
      
      String fileNameSolved = "batch/solved_generic_distribution_instances_MILP_receding.json";
      SKPGenericDistributionRecedingSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, linearizationSamples, simulationRuns);
      
      //solvedBatch = retrieveSolvedBatchMILP(fileNameSolved); // Batch cannot be retrieved because Distribution[] is not Serializable
      System.out.println(GSONUtility.<SKPGenericDistributionRecedingSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "batch/solved_generic_distribution_instances_MILP_receding.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPGenericDistributionRecedingSolvedInstance[] solveBatchMILP(SKPGenericDistribution[] instances, String fileName, int linearizationSamples, int simulationRuns) throws IloException {
      ArrayList<SKPGenericDistributionRecedingSolvedInstance>solved = new ArrayList<SKPGenericDistributionRecedingSolvedInstance>();
      for(SKPGenericDistribution instance : instances) {
         solved.add(new SimulateGenericDistributionReceding(instance).solve(simulationRuns, linearizationSamples));
         GSONUtility.<SKPGenericDistributionRecedingSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPGenericDistributionRecedingSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPGenericDistributionRecedingSolvedInstance[solved.size()]);
   }
   
   private static void storeSolvedBatchToCSV(SKPGenericDistributionRecedingSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights "
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
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}
