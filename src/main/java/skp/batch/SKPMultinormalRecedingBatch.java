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

import ilog.concert.IloException;

import skp.instance.SKPMultinormal;
import skp.sim.SimulateMultinormalReceding;
import skp.sim.instance.SKPMultinormalRecedingSolvedInstance;
import skp.utilities.gson.GSONUtility;

public class SKPMultinormalRecedingBatch extends SKPMultinormalBatch{
   
   public static void main(String args[]) {
      int[] instanceSize = {25, 50, 100, 500};
      double[] coeff_of_var  = {0.1, 0.2};
      double[] coeff_of_cor  = {0.75, 0.95};
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
               for(double rho : coeff_of_cor) {
                  for(boolean ignoreCorrelation : new boolean[] {false, true}) { // switch to ignore correlation while solving MILP model
                  
                     String batchFileName = "batch/"+t.toString()+"/"+size+"/"+cv+"/"+rho+"/multinormal_instances.json";
                     
                     /**
                      *  Generate instances using SKPMultinormalBatch
                      *  
                      *  generateInstances(batchFileName);
                      */
                     
                     int partitions = 10;
                     int simulationRuns = 100;
                     try {
                        solveMILP(batchFileName, partitions, simulationRuns, ignoreCorrelation, "batch/"+t.toString()+"/"+size+"/"+cv+"/"+rho);
                     } catch (IloException e) {
                        e.printStackTrace();
                     }
                  }
               }
            }
         }
      }
      
       
   }

   /*
    * MILP receding
    */ 
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns, boolean ignoreCorrelation, String folder) throws IloException {
      SKPMultinormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = ignoreCorrelation ? folder+"/solved_multinormal_instances_MILP_receding_ignore_correlation.json" : folder+"/solved_multinormal_instances_MILP_receding.json";
      SKPMultinormalRecedingSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns, ignoreCorrelation);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPMultinormalRecedingSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = ignoreCorrelation ? folder+"/solved_multinormal_instances_MILP_receding_ignore_correlation.csv" : folder+"/solved_multinormal_instances_MILP_receding.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPMultinormalRecedingSolvedInstance[] solveBatchMILP(SKPMultinormal[] instances, String fileName, int partitions, int simulationRuns, boolean ignoreCorrelation) throws IloException {
      ArrayList<SKPMultinormalRecedingSolvedInstance>solved = new ArrayList<SKPMultinormalRecedingSolvedInstance>();
      for(SKPMultinormal instance : instances) {
         solved.add(new SimulateMultinormalReceding(instance, partitions).solve(simulationRuns, ignoreCorrelation));
         System.out.println("Solved receding horizon instance number "+solved.size());
         GSONUtility.<SKPMultinormalRecedingSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPMultinormalRecedingSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPMultinormalRecedingSolvedInstance[solved.size()]);
   }
   
   private static void storeSolvedBatchToCSV(SKPMultinormalRecedingSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights, covarianceWeights, "
            + "capacity, shortageCost, simulatedSolutionMean, simulatedSolutionStd, "
            + "simulationRuns, piecewisePartitions, "
            + "piecewiseSamples, evp, evwpi,"
            + "evwpi_obj_2_n\n";
      String body = "";
      
      for(SKPMultinormalRecedingSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
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
         e.printStackTrace();
      }
   }
   
   private static SKPMultinormalRecedingSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPMultinormalRecedingSolvedInstance[] solvedInstances = GSONUtility.<SKPMultinormalRecedingSolvedInstance[]>retrieveJSONInstance(fileName, SKPMultinormalRecedingSolvedInstance[].class);
      return solvedInstances;
   }
}
