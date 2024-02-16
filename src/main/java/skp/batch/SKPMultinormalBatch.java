/**
 * To run from Mac OS
 * 
 * -Djava.library.path=/Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
 * 
 * Environment variable
 * 
 * DYLD_LIBRARY_PATH /Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
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
import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.SKPMultinormal;
import skp.milp.SKPMultinormalMILP;
import skp.milp.instance.SKPMultinormalMILPSolvedInstance;
import skp.utililities.gson.GSONUtility;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.probdist.UniformIntDist;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.randvar.RandomVariateGenInt;

public class SKPMultinormalBatch extends SKPBatch {
   
   public static void main(String args[]) {
      File folder = new File("scrap");
      if (!folder.exists()) {
        folder.mkdir();
      } 
      
      String batchFileName = "scrap/multinormal_instances.json";
      generateInstances(batchFileName);
      
      int partitions = 10;
      String OPLDataFileZipArchive = "scrap/multinormal_instances_opl.zip";
      storeBatchAsOPLDataFiles(retrieveBatch(batchFileName), OPLDataFileZipArchive, partitions);
      
      int simulationRuns = 100000;
      try {
         solveMILP(batchFileName, partitions, simulationRuns);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
   
   private static void generateInstances(String batchFileName) {
      int instances = 10;
      int instanceSize = 10;
      
      Distribution expectedValuePerUnit = new UniformDist(0.1,10);
      Distribution expectedWeight = new UniformDist(15,70);
      Distribution coefficientOfVariation = new UniformDist(0.1, 0.5);
      Distribution correlationCoefficient = new UniformDist(0, 1);
      DiscreteDistributionInt capacity = new UniformIntDist(100,200);
      Distribution shortageCost = new UniformDist(50,150);
      
      SKPMultinormalBatch batch = new SKPMultinormalBatch(expectedValuePerUnit, expectedWeight, coefficientOfVariation, correlationCoefficient, capacity, shortageCost);
      batch.generateBatch(instances, instanceSize, batchFileName);
   }
   
   protected Distribution coefficientOfVariation;
   protected Distribution correlationCoefficient;
   
   public SKPMultinormalBatch(
         Distribution expectedValuePerUnit,
         Distribution expectedWeight,
         Distribution coefficientOfVariation,
         Distribution correlationCoefficient,
         DiscreteDistributionInt capacity,
         Distribution shortageCost) {
      super(expectedValuePerUnit, expectedWeight, capacity, shortageCost);
      this.coefficientOfVariation = coefficientOfVariation;
      this.correlationCoefficient = correlationCoefficient;
   }
   
   /*
    * MILP
    */   
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns) throws IloException {
      SKPMultinormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "scrap/solvedMultinormalInstancesMILP.json";
      SKPMultinormalMILPSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPMultinormalMILPSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solvedMultinormalInstancesMILP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPMultinormalMILPSolvedInstance[] solveBatchMILP(SKPMultinormal[] instances, String fileName, int partitions, int simulationRuns) throws IloException {
      ArrayList<SKPMultinormalMILPSolvedInstance>solved = new ArrayList<SKPMultinormalMILPSolvedInstance>();
      for(SKPMultinormal instance : instances) {
         solved.add(new SKPMultinormalMILP(instance, partitions).solve(simulationRuns));
         GSONUtility.<SKPMultinormalMILPSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPMultinormalMILPSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPMultinormalMILPSolvedInstance[solved.size()]);
   }

   private static void storeSolvedBatchToCSV(SKPMultinormalMILPSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValuesPerUnit, expectedWeights, covarianceWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, milpOptimalityGap, piecewisePartitions, "
            + "piecewiseSamples, milpMaxLinearizationError, simulatedLinearizationError,"
            + "cplexSolutionTimeMs, simplexIterations, exploredNodes\n";
      String body = "";
      
      for(SKPMultinormalMILPSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(s.instance.getWeights().getMean()).replace(",", "\t")+ ", " +
                 Arrays.deepToString(s.instance.getWeights().getCovariance()).replace(",", "\t")+ ", " +
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
   
   private static SKPMultinormalMILPSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPMultinormalMILPSolvedInstance[] solvedInstances = GSONUtility.<SKPMultinormalMILPSolvedInstance[]>retrieveJSONInstance(fileName, SKPMultinormalMILPSolvedInstance[].class);
      return solvedInstances;
   }
   
   /**
    * Generate a batch of instances
    */
   
   public void generateBatch(int numberOfInstances, int instanceSize, String fileName) {
      SKPMultinormal[] instances = this.generateInstances(numberOfInstances, instanceSize);
      GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(instances, fileName);
   }
   
   public static SKPMultinormal[] retrieveBatch(String fileName) {
      SKPMultinormal[] instances = GSONUtility.<SKPMultinormal[]>retrieveJSONInstance(fileName, SKPMultinormal[].class);
      return instances;
   }
   
   private SKPMultinormal[] generateInstances(int numberOfInstances, int instanceSize){
      randGenerator.setSeed(seed);
      randGenerator.resetStartStream();
      SKPMultinormal[] instances = IntStream.iterate(0, i -> i + 1)
                                            .limit(numberOfInstances)
                                            .mapToObj(i -> new SKPMultinormal(
                                                  (new RandomVariateGen(randGenerator, this.expectedValuePerUnit)).nextArrayOfDouble(instanceSize),
                                                  (new RandomVariateGen(randGenerator, this.expectedWeight)).nextArrayOfDouble(instanceSize),
                                                  (new RandomVariateGen(randGenerator, this.coefficientOfVariation)).nextDouble(),
                                                  (new RandomVariateGen(randGenerator, this.correlationCoefficient)).nextDouble(),
                                                  (new RandomVariateGenInt(randGenerator, this.capacity)).nextInt(),
                                                  (new RandomVariateGen(randGenerator, this.shortageCost)).nextDouble()))
                                            .toArray(SKPMultinormal[]::new);
      return instances;
   }
   
   public static void storeBatchAsOPLDataFiles(SKPMultinormal[] instances, String OPLDataFileZipArchive, int partitions) {
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

         for(SKPMultinormal s : instances) {
            String body = "";
            body += 
                  "N = "+s.getItems()+";\n"+
                        "expectedValues = "+Arrays.toString(s.getExpectedValues())+";\n"+
                        "expectedWeights = "+Arrays.toString(s.getWeights().getMean())+";\n"+
                        "varianceCovarianceWeights = "+Arrays.deepToString(s.getWeights().getCovariance())+";\n"+
                        "C = "+s.getCapacity()+";\n"+
                        "c = "+s.getShortageCost()+";\n\n"+
                        "nbpartitions = "+partitions+";\n"+
                        "prob = "+Arrays.toString(PiecewiseStandardNormalFirstOrderLossFunction.getProbabilities(partitions))+";\n"+
                        "means = "+Arrays.toString(PiecewiseStandardNormalFirstOrderLossFunction.getMeans(partitions))+";\n"+
                        "error = "+PiecewiseStandardNormalFirstOrderLossFunction.getError(partitions)+";";

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
