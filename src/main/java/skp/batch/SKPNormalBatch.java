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
import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.SKPMultinormal;
import skp.instance.SKPNormal;
import skp.milp.SKPNormalMILP;
import skp.milp.instance.SKPNormalMILPSolvedInstance;
import skp.sdp.DSKPNormal;
import skp.sdp.instance.DSKPNormalSolvedInstance;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.probdist.UniformIntDist;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.randvar.RandomVariateGenInt;

public class SKPNormalBatch extends SKPBatch {

   public static void main(String args[]) {
      File folder = new File("batch");
      if (!folder.exists()) {
        folder.mkdir();
      } 
      
      String batchFileName = "scrap/normal_instances.json";
      String multinormalBatchFileName = "scrap/multinormal_instances.json";
      
      generateInstances(batchFileName, INSTANCE_TYPE.NORMAL);
      generateInstances(multinormalBatchFileName, INSTANCE_TYPE.MULTINORMAL);
      
      int partitions = 10;
      String OPLDataFileZipArchive = "scrap/normal_instances_opl.zip";
      storeBatchAsOPLDataFiles(retrieveBatch(batchFileName), OPLDataFileZipArchive, 10);
      
      int simulationRuns = 100000;
      try {
         solveMILP(batchFileName, partitions, simulationRuns);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      solveDSKP(batchFileName);
   }
   
   enum INSTANCE_TYPE {
      NORMAL,
      MULTINORMAL
   };
   
   private static void generateInstances(String batchFileName, INSTANCE_TYPE type) {
      int instances = 10;
      int instanceSize = 10;
      
      Distribution expectedValuePerUnit = new UniformDist(0.1,10);
      Distribution expectedWeight = new UniformDist(15,70);
      Distribution coefficientOfVariation = new UniformDist(0.1, 0.5);
      DiscreteDistributionInt capacity = new UniformIntDist(100,200);
      Distribution shortageCost = new UniformDist(50,150);
      
      SKPNormalBatch batch = new SKPNormalBatch(expectedValuePerUnit, expectedWeight, coefficientOfVariation, capacity, shortageCost);
      switch(type) {
         case NORMAL:
            batch.generateBatch(instances, instanceSize, batchFileName);
            break;
         case MULTINORMAL:
            batch.generateMultinormalBatch(instances, instanceSize, batchFileName);
            break;
      }
   }
   
   protected Distribution coefficientOfVariation;
   
   public SKPNormalBatch(
         Distribution expectedValuePerUnit,
         Distribution expectedWeight,
         Distribution coefficientOfVariation,
         DiscreteDistributionInt capacity,
         Distribution shortageCost) {
      super(expectedValuePerUnit, expectedWeight, capacity, shortageCost);
      this.coefficientOfVariation = coefficientOfVariation;
   }
   
   /*
    * MILP
    */   
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns) throws IloException {
      SKPNormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "scrap/solved_normal_instances_MILP.json";
      SKPNormalMILPSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPNormalMILPSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solved_normal_instances_MILP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static SKPNormalMILPSolvedInstance[] solveBatchMILP(SKPNormal[] instances, String fileName, int partitions, int simulationRuns) throws IloException {
      ArrayList<SKPNormalMILPSolvedInstance>solved = new ArrayList<SKPNormalMILPSolvedInstance>();
      for(SKPNormal instance : instances) {
         solved.add(new SKPNormalMILP(instance, partitions).solve(simulationRuns));
         GSONUtility.<SKPNormalMILPSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPNormalMILPSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPNormalMILPSolvedInstance[solved.size()]);
   }

   private static void storeSolvedBatchToCSV(SKPNormalMILPSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValuesPerUnit, expectedWeights, stdWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, milpOptimalityGap, piecewisePartitions, "
            + "piecewiseSamples, milpMaxLinearizationError, simulatedLinearizationError,"
            + "cplexSolutionTimeMs, simplexIterations, exploredNodes\n";
      String body = "";
      
      for(SKPNormalMILPSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).mapToDouble(d -> d.getMean()).toArray()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).mapToDouble(d -> d.getSigma()).toArray()).replace(",", "\t")+ ", " +
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
   
   private static SKPNormalMILPSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPNormalMILPSolvedInstance[] solvedInstances = GSONUtility.<SKPNormalMILPSolvedInstance[]>retrieveJSONInstance(fileName, SKPNormalMILPSolvedInstance[].class);
      return solvedInstances;
   }
   
   /*
    * DSKP
    */
   
   public static void solveDSKP(String fileName) {
      SKPNormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "scrap/solved_normal_instances_DSKP.json";
      DSKPNormalSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved);
      System.out.println(GSONUtility.<DSKPNormalSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solved_normal_instances_DSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static DSKPNormalSolvedInstance[] solveBatchDSKP(SKPNormal[] instances, String fileName) {
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPNormalSolvedInstance>solved = new ArrayList<DSKPNormalSolvedInstance>();
      for(SKPNormal instance : instances) {
         solved.add(new DSKPNormal(instance, truncationQuantile).solve());
         GSONUtility.<DSKPNormalSolvedInstance[]>saveInstanceToJSON(solved.toArray(new DSKPNormalSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPNormalSolvedInstance[solved.size()]);
   }

   private static void storeSolvedBatchToCSV(DSKPNormalSolvedInstance[] instances, String fileName) {
      String header = "instanceID, expectedValuesPerUnit, expectedWeights, stdWeights, capacity, shortageCost, solutionValue, solutionTimeMs, statesExplored\n";
      String body = "";
      
      for(DSKPNormalSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).mapToDouble(d -> d.getMean()).toArray()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).mapToDouble(d -> d.getSigma()).toArray()).replace(",", "\t")+ ", " +
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
   
   private static DSKPNormalSolvedInstance[] retrieveSolvedBatchDSKP(String fileName) {
      DSKPNormalSolvedInstance[] solvedInstances = GSONUtility.<DSKPNormalSolvedInstance[]>retrieveJSONInstance(fileName, DSKPNormalSolvedInstance[].class);
      return solvedInstances;
   }   
   
   /**
    * Generate a batch of instances
    */
   
   public void generateBatch(int numberOfInstances, int instanceSize, String fileName) {
      SKPNormal[] instances = this.generateInstances(numberOfInstances, instanceSize);
      GSONUtility.<SKPNormal[]>saveInstanceToJSON(instances, fileName);
   }
   
   public void generateMultinormalBatch(int numberOfInstances, int instanceSize, String fileName) {
      SKPMultinormal[] instances = this.generateMultinormalInstances(numberOfInstances, instanceSize);
      GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(instances, fileName);
   }
   
   public static SKPNormal[] retrieveBatch(String fileName) {
      SKPNormal[] instances = GSONUtility.<SKPNormal[]>retrieveJSONInstance(fileName, SKPNormal[].class);
      return instances;
   }
   
   private SKPNormal[] generateInstances(int numberOfInstances, int instanceSize){
      randGenerator.setSeed(seed);
      randGenerator.resetStartStream();
      SKPNormal[] instances = IntStream.iterate(0, i -> i + 1)
                                        .limit(numberOfInstances)
                                        .mapToObj(i -> new SKPNormal(
                                              (new RandomVariateGen(randGenerator, this.expectedValuePerUnit)).nextArrayOfDouble(instanceSize),
                                              (new RandomVariateGen(randGenerator, this.expectedWeight)).nextArrayOfDouble(instanceSize),
                                              (new RandomVariateGen(randGenerator, this.coefficientOfVariation)).nextDouble(),
                                              (new RandomVariateGenInt(randGenerator, this.capacity)).nextInt(),
                                              (new RandomVariateGen(randGenerator, this.shortageCost)).nextDouble()))
                                        .toArray(SKPNormal[]::new);
      return instances;
   }
   
   private SKPMultinormal[] generateMultinormalInstances(int numberOfInstances, int instanceSize){
      randGenerator.setSeed(seed);
      randGenerator.resetStartStream();
      SKPMultinormal[] instances = IntStream.iterate(0, i -> i + 1)
                                        .limit(numberOfInstances)
                                        .mapToObj(i -> new SKPMultinormal(
                                              (new RandomVariateGen(randGenerator, this.expectedValuePerUnit)).nextArrayOfDouble(instanceSize),
                                              (new RandomVariateGen(randGenerator, this.expectedWeight)).nextArrayOfDouble(instanceSize),
                                              (new RandomVariateGen(randGenerator, this.coefficientOfVariation)).nextDouble(),
                                              (new RandomVariateGenInt(randGenerator, this.capacity)).nextInt(),
                                              (new RandomVariateGen(randGenerator, this.shortageCost)).nextDouble()))
                                        .toArray(SKPMultinormal[]::new);
      return instances;
   }

   public static void storeBatchAsOPLDataFiles(SKPNormal[] instances, String OPLDataFileZipArchive, int partitions) {
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

         for(SKPNormal s : instances) {
            String body = "";
            body += 
                  "N = "+s.getItems()+";\n"+
                        "expectedValues = "+Arrays.toString(s.getExpectedValues())+";\n"+
                        "expectedWeights = "+Arrays.toString(Arrays.stream(s.getWeights()).mapToDouble(d -> d.getMean()).toArray())+";\n"+
                        "varianceWeights = "+Arrays.toString(Arrays.stream(s.getWeights()).mapToDouble(d -> d.getVariance()).toArray())+";\n"+
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
