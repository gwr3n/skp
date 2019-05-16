package skp.batch;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import ilog.concert.IloException;

import skp.instance.SKPMultinormal;
import skp.sim.SimulateMultinormalReceding;
import skp.sim.instance.SKPMultinormalRecedingSolvedInstance;
import skp.utililities.gson.GSONUtility;

public class SKPMultinormalRecedingBatch extends SKPMultinormalBatch{
   public static void main(String args[]) {
      String batchFileName = "scrap/multinormal_instances.json";
      int instances = 10;
      int instanceSize = 10;
      generateBatch(instances, instanceSize, batchFileName);
      
      int partitions = 10;
      String OPLDataFileZipArchive = "scrap/multinormal_instances_opl.zip";
      storeBatchAsOPLDataFiles(retrieveBatch(batchFileName), OPLDataFileZipArchive, partitions);
      
      int simulationRuns = 100;
      try {
         solveMILP(batchFileName, partitions, simulationRuns);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }

   /*
    * MILP receding
    */ 
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns) throws IloException {
      SKPMultinormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "scrap/solvedMultinormalInstancesMILPReceding.json";
      SKPMultinormalRecedingSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPMultinormalRecedingSolvedInstance[]>printInstanceAsGSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solvedMultinormalInstancesMILPReceding.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPMultinormalRecedingSolvedInstance[] solveBatchMILP(SKPMultinormal[] instances, String fileName, int partitions, int simulationRuns) throws IloException {
      ArrayList<SKPMultinormalRecedingSolvedInstance>solved = new ArrayList<SKPMultinormalRecedingSolvedInstance>();
      for(SKPMultinormal instance : instances) {
         solved.add(new SimulateMultinormalReceding(instance, partitions).solve(simulationRuns));
         GSONUtility.<SKPMultinormalRecedingSolvedInstance[]>saveInstanceToGSON(solved.toArray(new SKPMultinormalRecedingSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPMultinormalRecedingSolvedInstance[solved.size()]);
   }
   
   private static void storeSolvedBatchToCSV(SKPMultinormalRecedingSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValuesPerUnit, expectedWeights, covarianceWeights, "
            + "capacity, shortageCost, simulatedSolutionMean, simulatedSolutionStd, "
            + "simulationRuns, piecewisePartitions, "
            + "piecewiseSamples, evp, evwpi,"
            + "evwpi_obj_2_n\n";
      String body = "";
      
      for(SKPMultinormalRecedingSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(s.instance.getWeights().getMean()).replace(",", "\t")+ ", " +
                 Arrays.deepToString(s.instance.getWeights().getCovariance()).replace(",", "\t")+ ", " +
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
   
   private static SKPMultinormalRecedingSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPMultinormalRecedingSolvedInstance[] solvedInstances = GSONUtility.<SKPMultinormalRecedingSolvedInstance[]>retrieveInstance(fileName, SKPMultinormalRecedingSolvedInstance[].class);
      return solvedInstances;
   }
}
