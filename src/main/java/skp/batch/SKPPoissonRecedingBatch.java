package skp.batch;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import ilog.concert.IloException;

import skp.instance.SKPPoisson;
import skp.sim.SimulatePoissonReceding;
import skp.sim.instance.SKPPoissonRecedingSolvedInstance;
import skp.utilities.gson.GSONUtility;

public class SKPPoissonRecedingBatch extends SKPPoissonBatch{
   
   public static void main(String args[]) {
      File folder = new File("batch");
      if (!folder.exists()) {
        folder.mkdir();
      } 
      
      String batchFileName = "batch/poisson_instances.json";
      
      /**
       *  Generate instances using SKPPoissonBatch
       *  
       *  generateInstances(batchFileName);
       */
      
      int partitions = 10;
      int linearizationSamples = 1000;
      int simulationRuns = 100;
      try {
         solveMILP(batchFileName, partitions, simulationRuns, linearizationSamples);
      } catch (IloException e) {
         e.printStackTrace();
      }
   }

   /*
    * MILP receding
    */ 
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns, int linearizationSamples) throws IloException {
      SKPPoisson[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "batch/solved_poisson_instances_MILP_receding.json";
      SKPPoissonRecedingSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns, linearizationSamples);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPPoissonRecedingSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "batch/solved_poisson_instances_MILP_receding.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPPoissonRecedingSolvedInstance[] solveBatchMILP(SKPPoisson[] instances, String fileName, int partitions, int simulationRuns, int linearizationSamples) throws IloException {
      ArrayList<SKPPoissonRecedingSolvedInstance>solved = new ArrayList<SKPPoissonRecedingSolvedInstance>();
      for(SKPPoisson instance : instances) {
         solved.add(new SimulatePoissonReceding(instance, partitions).solve(simulationRuns, linearizationSamples));
         GSONUtility.<SKPPoissonRecedingSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPPoissonRecedingSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPPoissonRecedingSolvedInstance[solved.size()]);
   }
   
   private static void storeSolvedBatchToCSV(SKPPoissonRecedingSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights "
            + "capacity, shortageCost, simulatedSolutionMean, simulatedSolutionStd, "
            + "simulationRuns, piecewisePartitions, "
            + "piecewiseSamples, evp, evwpi,"
            + "evwpi_obj_2_n\n";
      String body = "";
      
      for(SKPPoissonRecedingSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).mapToDouble(d -> d.getMean()).toArray()).replace(",", "\t")+ ", " +
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
         e.printStackTrace();
      }
   }
   
   private static SKPPoissonRecedingSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPPoissonRecedingSolvedInstance[] solvedInstances = GSONUtility.<SKPPoissonRecedingSolvedInstance[]>retrieveJSONInstance(fileName, SKPPoissonRecedingSolvedInstance[].class);
      return solvedInstances;
   }
}
