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

public class SKPNormalRecedingBatch extends SKPNormalBatch{
   
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
               String batchFileName = "batch/"+t.toString()+"/"+size+"/"+cv+"/normal_instances.json";
               
               /**
                *  Generate instances using SKPNormalBatch
                *  
                *  generateInstances(batchFileName);
                */
               
               int partitions = 10;
               int simulationRuns = 100;
               try {
                  solveMILP(batchFileName, partitions, simulationRuns, "batch/"+t.toString()+"/"+size+"/"+cv);
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
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns, String folder) throws IloException {
      SKPNormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = folder+"/solved_normal_instances_MILP_receding.json";
      SKPNormalRecedingSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPNormalRecedingSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = folder+"/solved_normal_instances_MILP_receding.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPNormalRecedingSolvedInstance[] solveBatchMILP(SKPNormal[] instances, String fileName, int partitions, int simulationRuns) throws IloException {
      ArrayList<SKPNormalRecedingSolvedInstance>solved = new ArrayList<SKPNormalRecedingSolvedInstance>();
      for(SKPNormal instance : instances) {
         solved.add(new SimulateNormalReceding(instance, partitions).solve(simulationRuns));
         System.out.println("Solved instance number "+solved.size());
         GSONUtility.<SKPNormalRecedingSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPNormalRecedingSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPNormalRecedingSolvedInstance[solved.size()]);
   }
   
   private static void storeSolvedBatchToCSV(SKPNormalRecedingSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights, stdWeights, "
            + "capacity, shortageCost, simulatedSolutionMean, simulatedSolutionStd, "
            + "simulationRuns, piecewisePartitions, "
            + "piecewiseSamples, evp, evwpi,"
            + "evwpi_obj_2_n\n";
      String body = "";
      
      for(SKPNormalRecedingSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
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
         e.printStackTrace();
      }
   }
   
   private static SKPNormalRecedingSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPNormalRecedingSolvedInstance[] solvedInstances = GSONUtility.<SKPNormalRecedingSolvedInstance[]>retrieveJSONInstance(fileName, SKPNormalRecedingSolvedInstance[].class);
      return solvedInstances;
   }
}
