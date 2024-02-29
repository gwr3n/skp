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

import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;

public class SKPMultinormalRecedingBatch extends SKPMultinormalBatch{
   
   public static void main(String args[]) {
      File folder = new File("batch");
      if (!folder.exists()) {
        folder.mkdir();
      } 
      
      String batchFileName = "batch/multinormal_instances.json";
      
      /**
       *  Generate instances using SKPMultinormalBatch
       *  
       *  generateInstances(batchFileName);
       */
      
      int partitions = 10;
      int simulationRuns = 100;
      try {
         solveMILP(batchFileName, partitions, simulationRuns);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      } 
   }
     
   public SKPMultinormalRecedingBatch(
         Distribution expectedValuePerUnit,
         Distribution expectedWeight,
         Distribution coefficientOfVariation,
         Distribution correlationCoefficient,
         DiscreteDistributionInt capacity,
         Distribution shortageCost) {
      super(expectedValuePerUnit, expectedWeight, coefficientOfVariation, correlationCoefficient, capacity, shortageCost);
   }

   /*
    * MILP receding
    */ 
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns) throws IloException {
      SKPMultinormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "batch/solved_multinormal_instances_MILP_receding.json";
      SKPMultinormalRecedingSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPMultinormalRecedingSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "batch/solved_multinormal_instances_MILP_receding.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPMultinormalRecedingSolvedInstance[] solveBatchMILP(SKPMultinormal[] instances, String fileName, int partitions, int simulationRuns) throws IloException {
      ArrayList<SKPMultinormalRecedingSolvedInstance>solved = new ArrayList<SKPMultinormalRecedingSolvedInstance>();
      for(SKPMultinormal instance : instances) {
         solved.add(new SimulateMultinormalReceding(instance, partitions).solve(simulationRuns));
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
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
   
   private static SKPMultinormalRecedingSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPMultinormalRecedingSolvedInstance[] solvedInstances = GSONUtility.<SKPMultinormalRecedingSolvedInstance[]>retrieveJSONInstance(fileName, SKPMultinormalRecedingSolvedInstance[].class);
      return solvedInstances;
   }
}
