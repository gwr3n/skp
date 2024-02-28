package skp.batch;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import ilog.concert.IloException;

import skp.instance.SKPNormal;
import skp.sim.SimulateNormalReceding;
import skp.sim.instance.SKPNormalRecedingSolvedInstance;
import skp.utilities.gson.GSONUtility;

import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;

public class SKPNormalRecedingBatch extends SKPNormalBatch{
   
   public static void main(String args[]) {
      File folder = new File("batch");
      if (!folder.exists()) {
        folder.mkdir();
      } 
      
      String batchFileName = "batch/normal_instances.json";
      
      /**
       *  Generate instances using SKPNormalBatch
       *  
       *  generateInstances(batchFileName);
       */
      
      int partitions = 10;
      String OPLDataFileZipArchive = "batch/normal_instances_opl.zip";
      storeBatchAsOPLDataFiles(retrieveBatch(batchFileName), OPLDataFileZipArchive, partitions);
      
      int simulationRuns = 100;
      try {
         solveMILP(batchFileName, partitions, simulationRuns);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      } 
   }
     
   public SKPNormalRecedingBatch(
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
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns) throws IloException {
      SKPNormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "batch/solved_normal_instances_MILP_receding.json";
      SKPNormalRecedingSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPNormalRecedingSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "batch/solved_multinormal_instances_MILP_receding.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPNormalRecedingSolvedInstance[] solveBatchMILP(SKPNormal[] instances, String fileName, int partitions, int simulationRuns) throws IloException {
      ArrayList<SKPNormalRecedingSolvedInstance>solved = new ArrayList<SKPNormalRecedingSolvedInstance>();
      for(SKPNormal instance : instances) {
         solved.add(new SimulateNormalReceding(instance, partitions).solve(simulationRuns));
         GSONUtility.<SKPNormalRecedingSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPNormalRecedingSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPNormalRecedingSolvedInstance[solved.size()]);
   }
   
   private static void storeSolvedBatchToCSV(SKPNormalRecedingSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValuesPerUnit, expectedWeights, stdWeights, "
            + "capacity, shortageCost, simulatedSolutionMean, simulatedSolutionStd, "
            + "simulationRuns, piecewisePartitions, "
            + "piecewiseSamples, evp, evwpi,"
            + "evwpi_obj_2_n\n";
      String body = "";
      
      for(SKPNormalRecedingSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).mapToDouble(d -> d.getMean()).toArray()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).mapToDouble(d -> d.getSigma()).toArray()).replace(",", "\t")+ ", " +
                 s.instance.getCapacity()+ ", " +
                 s.instance.getShortageCost()+ ", " +
                 s.simulatedSolutionMean + ", " +
                 s.simulatedSolutionStd + ", " +
                 s.simulationRuns + ", " +
                 s.piecewisePartitions + ", " +
                 s.piecewiseSamples + ", " +
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
   
   private static SKPNormalRecedingSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPNormalRecedingSolvedInstance[] solvedInstances = GSONUtility.<SKPNormalRecedingSolvedInstance[]>retrieveJSONInstance(fileName, SKPNormalRecedingSolvedInstance[].class);
      return solvedInstances;
   }
}
