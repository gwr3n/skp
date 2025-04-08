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
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import ilog.concert.IloException;

import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.SKPGenericDistribution;
import skp.instance.SKPMultinormal;
import skp.instance.SKPNormal;
import skp.milp.SKPNormalMILP;
import skp.milp.instance.SKPGenericDistributionCutsMVNSolvedInstance;
import skp.milp.instance.SKPNormalMILPSolvedInstance;
import skp.saa.instance.SKPGenericDistributionSAASolvedInstance;
import skp.sdp.DSKPNormal;
import skp.sdp.instance.DSKPNormalSolvedInstance;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.randvar.RandomVariateGen;

public class SKPNormalBatch extends SKPBatch {

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
               File folder = new File("batch/"+t.toString()+"/"+size+"/"+cv);
               if (!folder.exists()) {
                  folder.mkdirs();
               }
               
               String batchFileName = "batch/"+t.toString()+"/"+size+"/"+cv+"/normal_instances.json";
               generateInstances(batchFileName, t, size, cv);
               
               String OPLDataFileZipArchive = "batch/"+t.toString()+"/"+size+"/"+cv+"/normal_instances_opl.zip";
               storeBatchAsOPLDataFiles(retrieveBatch(batchFileName), OPLDataFileZipArchive, 10);
               
               int partitions = 10;
               int simulationRuns = 100000;
               
               int maxCuts = 1000;
               try {
                  solveMILP(batchFileName, partitions, simulationRuns, maxCuts, "batch/"+t.toString()+"/"+size+"/"+cv, METHOD.PWLA);
                  solveMILP(batchFileName, partitions, simulationRuns, maxCuts, "batch/"+t.toString()+"/"+size+"/"+cv, METHOD.DCG);
                  if(size == instanceSize[0]) 
                     solveMILP(batchFileName, partitions, simulationRuns, maxCuts, "batch/"+t.toString()+"/"+size+"/"+cv, METHOD.SAA);
               } catch (IloException e) {
                  e.printStackTrace();
               }
               if(size == instanceSize[0]) 
                  solveDSKP(batchFileName, "batch/"+t.toString()+"/"+size+"/"+cv);
            }
         }
      }
   }
   
   enum INSTANCE_TYPE {
      P05_UNCORRELATED,
      P05_WEAKLY_CORRELATED,
      P05_STRONGLY_CORRELATED,
      P05_INVERSE_STRONGLY_CORRELATED,
      P05_ALMOST_STRONGLY_CORRELATED,
      P05_SUBSET_SUM,
      P05_UNCORRELATED_SIMILAR_WEIGHTS,
      P05_PROFIT_CEILING,
      P05_CIRCLE_INSTANCES
   };
   
   /**
    * Generate a batch of instances
    */
   
   private static void generateInstances(String batchFileName, INSTANCE_TYPE type, int instanceSize, double cv) {
      int H = 10;
      int R = 100;
      double shortageCost = 10;
      
      switch(type) {
         /*case NORMAL: {
            int instances = 10;
            int instanceSize = 10;
            
            Distribution expectedValue = new UniformDist(2.75,275);
            Distribution expectedWeight = new UniformDist(15,70);
            Distribution coefficientOfVariation = new UniformDist(0.1, 0.3);
            DiscreteDistributionInt capacity = new UniformIntDist(100,200);
            Distribution shortageCost = new UniformDist(2,10);
            
            SKPNormal[] batch = new SKPNormal[instances];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            batch = IntStream.iterate(0, i -> i + 1)
                             .limit(instances)
                             .mapToObj(i -> new SKPNormal(
                                                    (new RandomVariateGen(randGenerator, expectedValue)).nextArrayOfDouble(instanceSize),
                                                    (new RandomVariateGen(randGenerator, expectedWeight)).nextArrayOfDouble(instanceSize),
                                                    (new RandomVariateGen(randGenerator, coefficientOfVariation)).nextDouble(),
                                                    (new RandomVariateGenInt(randGenerator, capacity)).nextInt(),
                                                    (new RandomVariateGen(randGenerator, shortageCost)).nextDouble()))
                             .toArray(SKPNormal[]::new);
            
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }*/
         case P05_UNCORRELATED: {
            
            SKPNormal[] batch = new SKPNormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPNormal(
                     expectedValues,
                     expectedWeights,
                     cv,
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_WEAKLY_CORRELATED: {
            
            SKPNormal[] batch = new SKPNormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double U = Math.max(1.0, Arrays.stream(expectedWeights).map(v -> v - R/10.0).max().getAsDouble());
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(U,U+2*R/10.0)).nextArrayOfDouble(instanceSize);
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPNormal(
                     expectedValues,
                     expectedWeights,
                     cv,
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_STRONGLY_CORRELATED: {
            
            SKPNormal[] batch = new SKPNormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> v + R/10.0).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPNormal(
                     expectedValues,
                     expectedWeights,
                     cv,
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_INVERSE_STRONGLY_CORRELATED: {
            
            SKPNormal[] batch = new SKPNormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedWeights = Arrays.stream(expectedValues).map(v -> v + R/10.0).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPNormal(
                     expectedValues,
                     expectedWeights,
                     cv,
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_ALMOST_STRONGLY_CORRELATED: {
            
            SKPNormal[] batch = new SKPNormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> new RandomVariateGen(randGenerator, new UniformDist(v + R/10.0 - R/500.0, v + R/10.0 + R/500.0)).nextDouble()).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPNormal(
                     expectedValues,
                     expectedWeights,
                     cv,
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_SUBSET_SUM: {
            
            SKPNormal[] batch = new SKPNormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.copyOf(expectedWeights, instanceSize);
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPNormal(
                     expectedValues,
                     expectedWeights,
                     cv,
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_UNCORRELATED_SIMILAR_WEIGHTS: {
            
            SKPNormal[] batch = new SKPNormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(R,R+10)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPNormal(
                     expectedValues,
                     expectedWeights,
                     cv,
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_PROFIT_CEILING: {
            
            SKPNormal[] batch = new SKPNormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double d = 3;
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> d*Math.ceil(v/d)).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPNormal(
                     expectedValues,
                     expectedWeights,
                     cv,
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_CIRCLE_INSTANCES: {
            
            SKPNormal[] batch = new SKPNormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> 2*Math.sqrt(4*R*R - Math.pow(v - 2*R,2))/3).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPNormal(
                     expectedValues,
                     expectedWeights,
                     cv,
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPNormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
      }
   }
   
   public static SKPNormal[] retrieveBatch(String fileName) {
      SKPNormal[] instances = GSONUtility.<SKPNormal[]>retrieveJSONInstance(fileName, SKPNormal[].class);
      return instances;
   }
   
   public static SKPGenericDistribution[] convertToGenericDistributionBatch(SKPNormal[] normalBatch) {
      SKPGenericDistribution[] batch = new SKPGenericDistribution[normalBatch.length];
      for(int i = 0; i < normalBatch.length; i++) {
         batch[i] = new SKPGenericDistribution(normalBatch[i]);
      }
      return batch;
   }
   
   public static SKPMultinormal[] convertToMVNDistributionBatch(SKPNormal[] normalBatch) {
      SKPMultinormal[] batch = new SKPMultinormal[normalBatch.length];
      for(int i = 0; i < normalBatch.length; i++) {
         batch[i] = new SKPMultinormal(normalBatch[i]);
      }
      return batch;
   }
   
   enum METHOD {
      PWLA,
      DCG,
      SAA
   }
   
   /*
    * MILP
    */   
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns, int maxCuts, String folder, METHOD method) throws IloException {
      switch(method){
      case DCG: //Compute optimal solution using Dynamic Cut Generation
         /*{
            SKPGenericDistribution[] batch = convertToGenericDistributionBatch(retrieveBatch(fileName));
            
            String fileNameSolved = folder+"/solved_normal_instances_DCG.json";
            SKPGenericDistributionCutsSolvedInstance[] solvedBatch = SKPGenericDistributionBatch.solveBatchMILPIterativeCuts(batch, fileNameSolved, linearizationSamples, maxCuts, simulationRuns);
            
            System.out.println(GSONUtility.<SKPGenericDistributionCutsSolvedInstance[]>printInstanceAsJSON(solvedBatch));
            
            String fileNameSolvedCSV = folder+"/solved_normal_instances_DCG.csv";
            SKPGenericDistributionBatch.storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
         }*/
         {
            SKPMultinormal[] batch = convertToMVNDistributionBatch(retrieveBatch(fileName));
            
            String fileNameSolved = folder+"/solved_normal_instances_DCG.json";
            SKPGenericDistributionCutsMVNSolvedInstance[] solvedBatch = SKPMultinormalBatch.solveBatchMILPDynamicCutGeneration(batch, fileNameSolved, maxCuts, simulationRuns);
            
            System.out.println(GSONUtility.<SKPGenericDistributionCutsMVNSolvedInstance[]>printInstanceAsJSON(solvedBatch));
            
            String fileNameSolvedCSV = folder+"/solved_normal_instances_DCG.csv";
            SKPGenericDistributionBatch.storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
         }
         break;
      case SAA: // SAA
         {  
            int Nsmall = 1000; 
            int Nlarge = simulationRuns; 
            int M = 1000;
            
            SKPGenericDistribution[] batch = convertToGenericDistributionBatch(retrieveBatch(fileName));
            
            String fileNameSolved = folder+"/solved_normal_instances_SAA.json";
            SKPGenericDistributionSAASolvedInstance[] solvedBatch = SKPGenericDistributionSAABatch.solveBatchMILPIterativeCuts(batch, fileNameSolved, Nsmall, Nlarge, M);
            
            System.out.println(GSONUtility.<SKPGenericDistributionSAASolvedInstance[]>printInstanceAsJSON(solvedBatch));
            
            String fileNameSolvedCSV = folder+"/solved_normal_instances_SAA.csv";
            SKPGenericDistributionSAABatch.storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
         }
         break;
      case PWLA:
      default:
         {
            SKPNormal[] batch = retrieveBatch(fileName);
            
            String fileNameSolved = folder+"/solved_normal_instances_MILP.json";
            SKPNormalMILPSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns);
            
            solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
            System.out.println(GSONUtility.<SKPNormalMILPSolvedInstance[]>printInstanceAsJSON(solvedBatch));
            
            String fileNameSolvedCSV = folder+"/solved_normal_instances_MILP.csv";
            storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
         }
      }
   }
   
   private static SKPNormalMILPSolvedInstance[] solveBatchMILP(SKPNormal[] instances, String fileName, int partitions, int simulationRuns) throws IloException {
      /*
       * Sequential
       *
      ArrayList<SKPNormalMILPSolvedInstance>solved = new ArrayList<SKPNormalMILPSolvedInstance>();
      for(SKPNormal instance : instances) {
         solved.add(new SKPNormalMILP(instance, partitions).solve(simulationRuns));
         GSONUtility.<SKPNormalMILPSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPNormalMILPSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPNormalMILPSolvedInstance[solved.size()]);*/
      
      /*
       * Parallel
       */
      SKPNormalMILPSolvedInstance[] solved = Arrays.stream(instances)
                                                   .parallel()
                                                   .map(instance -> {
                                                      try {
                                                         return SKPNormalMILP.solve(instance, partitions, simulationRuns);
                                                      } catch (IloException e) {
                                                         // TODO Auto-generated catch block
                                                         e.printStackTrace();
                                                         return null;
                                                      }
                                                   })
                                                   .toArray(SKPNormalMILPSolvedInstance[]::new);
      GSONUtility.<SKPNormalMILPSolvedInstance[]>saveInstanceToJSON(solved, fileName);
      return solved;
   }

   private static void storeSolvedBatchToCSV(SKPNormalMILPSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights, stdWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, milpOptimalityGap, piecewisePartitions, "
            + "piecewiseSamples, milpMaxLinearizationError, simulatedLinearizationError,"
            + "cplexSolutionTimeMs, simplexIterations, exploredNodes\n";
      String body = "";
      
      for(SKPNormalMILPSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
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
   
   public static void solveDSKP(String fileName, String folder) {
      SKPNormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = folder+"/solved_normal_instances_DSKP.json";
      DSKPNormalSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved);
      System.out.println(GSONUtility.<DSKPNormalSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = folder+"/solved_normal_instances_DSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static DSKPNormalSolvedInstance[] solveBatchDSKP(SKPNormal[] instances, String fileName) {
      /*
       * Sequential
       *
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPNormalSolvedInstance>solved = new ArrayList<DSKPNormalSolvedInstance>();
      for(SKPNormal instance : instances) {
         solved.add(new DSKPNormal(instance, truncationQuantile).solve());
         System.out.println("Solved DSKP instance number "+solved.size());
         GSONUtility.<DSKPNormalSolvedInstance[]>saveInstanceToJSON(solved.toArray(new DSKPNormalSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPNormalSolvedInstance[solved.size()]);*/
      
      /*
       * Parallel
       */
      double truncationQuantile = 0.999999999999999;
      DSKPNormalSolvedInstance[] solved = Arrays.stream(instances)
                                                .parallel()
                                                .map(instance -> new DSKPNormal(instance, truncationQuantile).solve())
                                                .toArray(DSKPNormalSolvedInstance[]::new);
      GSONUtility.<DSKPNormalSolvedInstance[]>saveInstanceToJSON(solved, fileName);
      return solved;
   }

   private static void storeSolvedBatchToCSV(DSKPNormalSolvedInstance[] instances, String fileName) {
      String header = "instanceID, expectedValues, expectedWeights, stdWeights, capacity, shortageCost, solutionValue, solutionTimeMs, statesExplored\n";
      String body = "";
      
      for(DSKPNormalSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
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
         e.printStackTrace();
      }
   }
   
   private static DSKPNormalSolvedInstance[] retrieveSolvedBatchDSKP(String fileName) {
      DSKPNormalSolvedInstance[] solvedInstances = GSONUtility.<DSKPNormalSolvedInstance[]>retrieveJSONInstance(fileName, DSKPNormalSolvedInstance[].class);
      return solvedInstances;
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
