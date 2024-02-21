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
import skp.sdp.DSKPMultinormal;
import skp.sdp.instance.DSKPMultinormalSolvedInstance;
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
       *  Generate instances using SKPMultinormalBatch or SKPNormalBatch
       *  
       *  generateInstances(batchFileName);
       */
      
      int partitions = 10;
      String OPLDataFileZipArchive = "batch/multinormal_instances_opl.zip";
      storeBatchAsOPLDataFiles(retrieveBatch(batchFileName), OPLDataFileZipArchive, partitions);
      
      int simulationRuns = 100;
      try {
         solveMILP(batchFileName, partitions, simulationRuns);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      /*
       * Note that this will only work if an instance covariance matrix takes the special structure 
       * $\rho^{|j-i|}\sigma_i\sigma_j$ discussed in [1], which ensures $P(d_t=x|d_{t-1}=y) = P(d_t=x|d_{t-1}=y,d_{t-2}=z,...)$.
       *
       * [1] M. Xiang, R. Rossi, B. Martin-Barragan, S. A. Tarim, "<a href="https://doi.org/10.1016/j.ejor.2022.04.011">
       * A mathematical programming-based solution method for the nonstationary inventory problem under correlated demand</a>," 
       * European Journal of Operational Research, Elsevier, Vol. 304(2): 515–524, 2023 
       */
      solveDSKP(batchFileName); 
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
      SKPMultinormalRecedingSolvedInstance[] solvedInstances = GSONUtility.<SKPMultinormalRecedingSolvedInstance[]>retrieveJSONInstance(fileName, SKPMultinormalRecedingSolvedInstance[].class);
      return solvedInstances;
   }
   
   /*
    * DSKP
    */
   
   public static void solveDSKP(String fileName) {
      SKPMultinormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "batch/solved_multinormal_instances_DSKP.json";
      DSKPMultinormalSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved);
      System.out.println(GSONUtility.<DSKPMultinormalSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "batch/solved_multinormal_instances_DSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static DSKPMultinormalSolvedInstance[] solveBatchDSKP(SKPMultinormal[] instances, String fileName) {
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPMultinormalSolvedInstance>solved = new ArrayList<DSKPMultinormalSolvedInstance>();
      for(SKPMultinormal instance : instances) {
         solved.add(new DSKPMultinormal(instance, truncationQuantile).solve());
         GSONUtility.<DSKPMultinormalSolvedInstance[]>saveInstanceToJSON(solved.toArray(new DSKPMultinormalSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPMultinormalSolvedInstance[solved.size()]);
   }

   private static void storeSolvedBatchToCSV(DSKPMultinormalSolvedInstance[] instances, String fileName) {
      String header = "instanceID, expectedValuesPerUnit, expectedWeights, covarianceWeights, capacity, shortageCost, solutionValue, solutionTimeMs, statesExplored\n";
      String body = "";
      
      for(DSKPMultinormalSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights().getMean()).toArray()).replace(",", "\t")+ ", " +
                 Arrays.deepToString(s.instance.getWeights().getCovariance()).replace(",", "\t")+ ", " +
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
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
   
   private static DSKPMultinormalSolvedInstance[] retrieveSolvedBatchDSKP(String fileName) {
      DSKPMultinormalSolvedInstance[] solvedInstances = GSONUtility.<DSKPMultinormalSolvedInstance[]>retrieveJSONInstance(fileName, DSKPMultinormalSolvedInstance[].class);
      return solvedInstances;
   }
}
