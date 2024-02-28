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
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.stream.IntStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import ilog.concert.IloException;
import skp.folf.PiecewiseFirstOrderLossFunction;
import skp.instance.SKPPoisson;
import skp.milp.SKPPoissonMILP;
import skp.milp.instance.SKPPoissonMILPSolvedInstance;
import skp.sdp.DSKPPoisson;
import skp.sdp.instance.DSKPPoissonSolvedInstance;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.probdist.UniformIntDist;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.randvar.RandomVariateGenInt;

public class SKPPoissonBatch extends SKPBatch {
   
   public static void main(String args[]) {
      File folder = new File("batch");
      if (!folder.exists()) {
        folder.mkdir();
      } 
      
      String batchFileName = "batch/poisson_instances.json";
      generateInstances(batchFileName);
      
      String OPLDataFileZipArchive = "batch/poisson_instances_opl.zip";
      int partitions = 10;
      int linearizationSamples = 1000;
      storeBatchAsOPLDataFiles(retrieveBatch(batchFileName), OPLDataFileZipArchive, partitions, linearizationSamples);
      
      int simulationRuns = 100000;
      try {
         solveMILP(batchFileName, partitions, linearizationSamples, simulationRuns);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      solveDSKP(batchFileName);
   }
   
   private static void generateInstances(String batchFileName) {
      int instances = 10;
      int instanceSize = 10;
      
      Distribution expectedValuePerUnit = new UniformDist(0.1,10);
      Distribution expectedWeight = new UniformDist(15,70);
      DiscreteDistributionInt capacity = new UniformIntDist(100,200);
      Distribution shortageCost = new UniformDist(50,150);
      
      SKPPoissonBatch batch = new SKPPoissonBatch(expectedValuePerUnit, expectedWeight, capacity, shortageCost);
      batch.generateBatch(instances, instanceSize, batchFileName);
   }
   
   public SKPPoissonBatch(
         Distribution expectedValuePerUnit,
         Distribution expectedWeight,
         DiscreteDistributionInt capacity,
         Distribution shortageCost) {
      super(expectedValuePerUnit, expectedWeight, capacity, shortageCost);
   }
   
   /*
    * MILP
    */
   
   public static void solveMILP(String fileName, int partitions, int linearizationSamples, int simulationRuns) throws IloException {
      SKPPoisson[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "batch/solved_poisson_instances_MILP.json";
      SKPPoissonMILPSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, linearizationSamples, simulationRuns);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPPoissonMILPSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "batch/solved_poisson_instances_MILP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPPoissonMILPSolvedInstance[] solveBatchMILP(SKPPoisson[] instances, String fileName, int partitions, int linearizationSamples, int simulationRuns) throws IloException {      
      ArrayList<SKPPoissonMILPSolvedInstance>solved = new ArrayList<SKPPoissonMILPSolvedInstance>();
      for(SKPPoisson instance : instances) {
         solved.add(new SKPPoissonMILP(instance, partitions, linearizationSamples).solve(simulationRuns));
         GSONUtility.<SKPPoissonMILPSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPPoissonMILPSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPPoissonMILPSolvedInstance[solved.size()]);
   }

   private static void storeSolvedBatchToCSV(SKPPoissonMILPSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValuesPerUnit, expectedWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, milpOptimalityGap, piecewisePartitions, "
            + "piecewiseSamples, milpMaxLinearizationError, simulatedLinearizationError,"
            + "cplexSolutionTimeMs, simplexIterations, exploredNodes\n";
      String body = "";
      
      for(SKPPoissonMILPSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).mapToDouble(d -> d.getMean()).toArray()).replace(",", "\t")+ ", " +
                 s.instance.getCapacity()+ ", " +
                 s.instance.getShortageCost()+ ", " +
                 Arrays.toString(s.optimalKnapsack).replace(",", "\t")+ ", " +
                 s.simulatedSolutionValue + ", " +
                 s.simulationRuns + ", " +
                 s.milpSolutionValue + ", " +
                 s.milpOptimalityGap + ", " +
                 s.piecewisePartitions + ", " +
                 s.piecewiseSamples + ", " +
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
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
   
   private static SKPPoissonMILPSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPPoissonMILPSolvedInstance[] solvedInstances = GSONUtility.<SKPPoissonMILPSolvedInstance[]>retrieveJSONInstance(fileName, SKPPoissonMILPSolvedInstance[].class);
      return solvedInstances;
   }
   
   /*
    * DSKP
    */
   
   public static void solveDSKP(String fileName) {
      SKPPoisson[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "batch/solved_poisson_instances_DSKP.json";
      DSKPPoissonSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved);
      System.out.println(GSONUtility.<DSKPPoissonSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "batch/solved_poisson_instances_DSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static DSKPPoissonSolvedInstance[] solveBatchDSKP(SKPPoisson[] instances, String fileName) {
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPPoissonSolvedInstance>solved = new ArrayList<DSKPPoissonSolvedInstance>();
      for(SKPPoisson instance : instances) {
         solved.add(new DSKPPoisson(instance, truncationQuantile).solve());
         GSONUtility.<DSKPPoissonSolvedInstance[]>saveInstanceToJSON(solved.toArray(new DSKPPoissonSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPPoissonSolvedInstance[solved.size()]);
   }

   private static void storeSolvedBatchToCSV(DSKPPoissonSolvedInstance[] instances, String fileName) {
      String header = "instanceID, expectedValuesPerUnit, expectedWeights, capacity, shortageCost, solutionValue, solutionTimeMs, statesExplored\n";
      String body = "";
      
      for(DSKPPoissonSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).mapToDouble(d -> d.getMean()).toArray()).replace(",", "\t")+ ", " +
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
   
   private static DSKPPoissonSolvedInstance[] retrieveSolvedBatchDSKP(String fileName) {
      DSKPPoissonSolvedInstance[] solvedInstances = GSONUtility.<DSKPPoissonSolvedInstance[]>retrieveJSONInstance(fileName, DSKPPoissonSolvedInstance[].class);
      return solvedInstances;
   }
   
   /**
    * Generate a batch of instances
    */
   
   public void generateBatch(int numberOfInstances, int instanceSize, String fileName) {
      SKPPoisson[] instances = this.generateInstances(numberOfInstances, instanceSize);
      GSONUtility.<SKPPoisson[]>saveInstanceToJSON(instances, fileName);
   }
   
   public static SKPPoisson[] retrieveBatch(String fileName) {
      SKPPoisson[] instances = GSONUtility.<SKPPoisson[]>retrieveJSONInstance(fileName, SKPPoisson[].class);
      return instances;
   }
   
   private SKPPoisson[] generateInstances(int numberOfInstances, int instanceSize){
      randGenerator.setSeed(seed);
      randGenerator.resetStartStream();
      SKPPoisson[] instances = IntStream.iterate(0, i -> i + 1)
                                        .limit(numberOfInstances)
                                        .mapToObj(i -> new SKPPoisson(
                                              (new RandomVariateGen(randGenerator, this.expectedValuePerUnit)).nextArrayOfDouble(instanceSize),
                                              (new RandomVariateGen(randGenerator, this.expectedWeight)).nextArrayOfDouble(instanceSize),
                                              (new RandomVariateGenInt(randGenerator, this.capacity)).nextInt(),
                                              (new RandomVariateGen(randGenerator, this.shortageCost)).nextDouble()))
                                        .toArray(SKPPoisson[]::new);
      return instances;
   }   
   
   protected static void storeBatchAsOPLDataFiles(SKPPoisson[] instances, String OPLDataFileZipArchive, int partitions, int linearizationSamples) {
      Date date = Calendar.getInstance().getTime();
      DateFormat dateFormat = new SimpleDateFormat("yyyy-mm-dd hh:mm:ss");
      String strDate = dateFormat.format(date);
      String header = 
            "/*********************************************\n" + 
            " * OPL 12.8.0.0 Data\n" + 
            " * Author: Roberto Rossi\n" + 
            " * Creation Date: "+strDate+"\n" + 
            " *********************************************/\n\n";
      
      try {
         File zipFile = new File(OPLDataFileZipArchive);
         ZipOutputStream out = new ZipOutputStream(new FileOutputStream(zipFile));

         for(SKPPoisson s : instances) {
            double[] probabilityMasses = new double[partitions];
            Arrays.fill(probabilityMasses, 1.0/partitions);
            
            String body = "";
            body += 
                  "N = "+s.getItems()+";\n"+
                        "expectedValues = "+Arrays.toString(s.getExpectedValues())+";\n"+
                        "expectedWeights = "+Arrays.toString(Arrays.stream(s.getWeights()).mapToDouble(d -> d.getMean()).toArray())+";\n"+
                        "C = "+s.getCapacity()+";\n"+
                        "c = "+s.getShortageCost()+";\n\n"+
                        "nbpartitions = "+partitions+";\n"+
                        "prob = "+Arrays.toString(probabilityMasses)+";\n"+
                        "means = "+Arrays.deepToString(PiecewiseFirstOrderLossFunction.poissonKnapsackPiecewiseFOLFConditionalExpectations(s.getCapacity(), probabilityMasses, linearizationSamples))+";\n"+
                        "error = "+Arrays.toString(PiecewiseFirstOrderLossFunction.poissonKnapsackPiecewiseFOLFApproximationErrors(s.getCapacity(), probabilityMasses, linearizationSamples))+";";

            ZipEntry e = new ZipEntry(s.getInstanceID()+".dat");
            out.putNextEntry(e);

            PrintWriter pw = new PrintWriter(out);
            pw.print(header+body);
            pw.flush();
            out.closeEntry();
         }
         out.close();
      }catch(IOException e) {
         e.printStackTrace();
      }
   }
}
