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
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import ilog.concert.IloException;
import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.SKPMultinormal;
import skp.milp.SKPMultinormalMILP;
import skp.milp.instance.SKPMultinormalMILPSolvedInstance;
import skp.sdp.DSKPMultinormal;
import skp.sdp.instance.DSKPMultinormalSolvedInstance;
import skp.utilities.gson.GSONUtility;

import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.randvar.RandomVariateGen;

public class SKPMultinormalBatch extends SKPBatch {
   
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
                  File folder = new File("batch/"+t.toString()+"/"+size+"/"+cv+"/"+rho);
                  if (!folder.exists()) {
                     folder.mkdirs();
                  }
                  
                  String batchFileName = "batch/"+t.toString()+"/"+size+"/"+cv+"/"+rho+"/multinormal_instances.json";
                  generateInstances(batchFileName, t, size, cv, rho);
                  
                  int partitions = 10;
                  String OPLDataFileZipArchive = "batch/"+t.toString()+"/"+size+"/"+cv+"/"+rho+"/multinormal_instances_opl.zip";
                  storeBatchAsOPLDataFiles(retrieveBatch(batchFileName), OPLDataFileZipArchive, partitions);
                  
                  int simulationRuns = 100000;
                  try {
                     solveMILP(batchFileName, partitions, simulationRuns, "batch/"+t.toString()+"/"+size+"/"+cv+"/"+rho);
                  } catch (IloException e) {
                     e.printStackTrace();
                  }
                  /*
                   * Note that solveDSKP will only work if an instance covariance matrix takes the special structure 
                   * $\rho^{|j-i|}\sigma_i\sigma_j$ discussed in [1], which ensures $P(d_t=x|d_{t-1}=y) = P(d_t=x|d_{t-1}=y,d_{t-2}=z,...)$.
                   *
                   * [1] M. Xiang, R. Rossi, B. Martin-Barragan, S. A. Tarim, "<a href="https://doi.org/10.1016/j.ejor.2022.04.011">
                   * A mathematical programming-based solution method for the nonstationary inventory problem under correlated demand</a>," 
                   * European Journal of Operational Research, Elsevier, Vol. 304(2): 515–524, 2023 
                   */
                  if(size == instanceSize[0])
                     solveDSKP(batchFileName, "batch/"+t.toString()+"/"+size+"/"+cv+"/"+rho);
               }
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
      
   private static void generateInstances(String batchFileName, INSTANCE_TYPE type, int instanceSize, double cv, double rho) {
      int H = 10;
      int R = 10;
      double shortageCost = 10;
      
      switch(type) {
         /*case MULTINORMAL: {
            int instances = 10;
                       
            Distribution expectedValue = new UniformDist(75,275);
            Distribution expectedWeight = new UniformDist(15,70);
            Distribution coefficientOfVariation = new UniformDist(0.1, 0.3);
            Distribution correlationCoefficient = new UniformDist(0, 1);
            DiscreteDistributionInt capacity = new UniformIntDist(100,200);
            Distribution shortageCost = new UniformDist(2,10);
            
            SKPMultinormal[] batch = new SKPMultinormal[instances];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            batch = IntStream.iterate(0, i -> i + 1)
                             .limit(instances)
                             .mapToObj(i -> new SKPMultinormal(
                                                        (new RandomVariateGen(randGenerator, expectedValue)).nextArrayOfDouble(instanceSize),
                                                        (new RandomVariateGen(randGenerator, expectedWeight)).nextArrayOfDouble(instanceSize),
                                                        (new RandomVariateGen(randGenerator, coefficientOfVariation)).nextDouble(),
                                                        (new RandomVariateGen(randGenerator, correlationCoefficient)).nextDouble(),
                                                        (new RandomVariateGenInt(randGenerator, capacity)).nextInt(),
                                                        (new RandomVariateGen(randGenerator, shortageCost)).nextDouble()))
                             .toArray(SKPMultinormal[]::new);
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }*/
         /**
          * Generate a batch of instances whose covariance matrix takes the special structure $\rho^{|j-i|}\sigma_i\sigma_j$ 
          * discussed in [1], which ensures $P(d_t=x|d_{t-1}=y) = P(d_t=x|d_{t-1}=y,d_{t-2}=z,...)$
          * 
          * [1] M. Xiang, R. Rossi, B. Martin-Barragan, S. A. Tarim, "<a href="https://doi.org/10.1016/j.ejor.2022.04.011">
          * A mathematical programming-based solution method for the nonstationary inventory problem under correlated demand</a>,
          * " European Journal of Operational Research, Elsevier, Vol. 304(2): 515–524, 2023
          */
         /*case SPECIAL_STRUCTURE: {
            int instances = 10;
            
            Distribution expectedValue = new UniformDist(2.75,275);
            Distribution expectedWeight = new UniformDist(5,10);
            Distribution coefficientOfVariation = new UniformDist(0.1, 0.3);
            Distribution correlationCoefficient = new UniformDist(0, 1);
            DiscreteDistributionInt capacity = new UniformIntDist(25,50);
            Distribution shortageCost = new UniformDist(2,10);
            
            SKPMultinormal[] batch = new SKPMultinormal[instances];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            batch = IntStream.iterate(0, i -> i + 1)
                             .limit(instances)
                             .mapToObj(i -> {double[] EV = (new RandomVariateGen(randGenerator, expectedValue)).nextArrayOfDouble(instanceSize);
                                             double[] EW = (new RandomVariateGen(randGenerator, expectedWeight)).nextArrayOfDouble(instanceSize);
                                             return new SKPMultinormal(
                                                     EV,
                                                     EW,
                                                     SKPMultinormal.calculateCovarianceSpecialStructure(EW, 
                                                                                         (new RandomVariateGen(randGenerator, coefficientOfVariation)).nextDouble(), 
                                                                                         (new RandomVariateGen(randGenerator, correlationCoefficient)).nextDouble()),
                                                     (new RandomVariateGenInt(randGenerator, capacity)).nextInt(),
                                                     (new RandomVariateGen(randGenerator, shortageCost)).nextDouble());})
                             .toArray(SKPMultinormal[]::new);
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }*/
         /**
          * Generate a batch of instances whose covariance matrix takes the special structure $\rho^{|j-i|}\sigma_i\sigma_j$ 
          * discussed in [1], which ensures $P(d_t=x|d_{t-1}=y) = P(d_t=x|d_{t-1}=y,d_{t-2}=z,...)$
          * 
          * [1] M. Xiang, R. Rossi, B. Martin-Barragan, S. A. Tarim, "<a href="https://doi.org/10.1016/j.ejor.2022.04.011">
          * A mathematical programming-based solution method for the nonstationary inventory problem under correlated demand</a>,
          * " European Journal of Operational Research, Elsevier, Vol. 304(2): 515–524, 2023
          */
         case P05_UNCORRELATED: {
            
            SKPMultinormal[] batch = new SKPMultinormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPMultinormal(
                     expectedValues,
                     expectedWeights,
                     SKPMultinormal.calculateCovarianceSpecialStructure(expectedWeights,cv,rho),
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_WEAKLY_CORRELATED: {
            
            SKPMultinormal[] batch = new SKPMultinormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double U = Math.max(1.0, Arrays.stream(expectedWeights).map(v -> v - R/10.0).max().getAsDouble());
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(U,U+2*R/10.0)).nextArrayOfDouble(instanceSize);
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPMultinormal(
                     expectedValues,
                     expectedWeights,
                     SKPMultinormal.calculateCovarianceSpecialStructure(expectedWeights,cv,rho),
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_STRONGLY_CORRELATED: {
            
            SKPMultinormal[] batch = new SKPMultinormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> v + R/10.0).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPMultinormal(
                     expectedValues,
                     expectedWeights,
                     SKPMultinormal.calculateCovarianceSpecialStructure(expectedWeights,cv,rho),
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_INVERSE_STRONGLY_CORRELATED: {
            
            SKPMultinormal[] batch = new SKPMultinormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedWeights = Arrays.stream(expectedValues).map(v -> v + R/10.0).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPMultinormal(
                     expectedValues,
                     expectedWeights,
                     SKPMultinormal.calculateCovarianceSpecialStructure(expectedWeights,cv,rho),
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_ALMOST_STRONGLY_CORRELATED: {
            
            SKPMultinormal[] batch = new SKPMultinormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> new RandomVariateGen(randGenerator, new UniformDist(v + R/10.0 - R/500.0, v + R/10.0 + R/500.0)).nextDouble()).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPMultinormal(
                     expectedValues,
                     expectedWeights,
                     SKPMultinormal.calculateCovarianceSpecialStructure(expectedWeights,cv,rho),
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_SUBSET_SUM: {
            
            SKPMultinormal[] batch = new SKPMultinormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.copyOf(expectedWeights, instanceSize);
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPMultinormal(
                     expectedValues,
                     expectedWeights,
                     SKPMultinormal.calculateCovarianceSpecialStructure(expectedWeights,cv,rho),
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_UNCORRELATED_SIMILAR_WEIGHTS: {
            
            SKPMultinormal[] batch = new SKPMultinormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(R,R+10)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPMultinormal(
                     expectedValues,
                     expectedWeights,
                     SKPMultinormal.calculateCovarianceSpecialStructure(expectedWeights,cv,rho),
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_PROFIT_CEILING: {
            
            SKPMultinormal[] batch = new SKPMultinormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double d = 3;
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> d*Math.ceil(v/d)).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPMultinormal(
                     expectedValues,
                     expectedWeights,
                     SKPMultinormal.calculateCovarianceSpecialStructure(expectedWeights,cv,rho),
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
         case P05_CIRCLE_INSTANCES: {
            
            SKPMultinormal[] batch = new SKPMultinormal[H];
            
            randGenerator.setSeed(seed);
            randGenerator.resetStartStream();
            
            for(int i = 0; i < H; i++) {
               double[] expectedWeights = new RandomVariateGen(randGenerator, new UniformDist(1,R)).nextArrayOfDouble(instanceSize);
               double[] expectedValues = Arrays.stream(expectedWeights).map(v -> 2*Math.sqrt(4*R*R - Math.pow(v - 2*R,2))/3).toArray();
               
               double capacity = ((i+1.0)/(H+1))*Arrays.stream(expectedWeights).sum();
               batch[i] = new SKPMultinormal(
                     expectedValues,
                     expectedWeights,
                     SKPMultinormal.calculateCovarianceSpecialStructure(expectedWeights,cv,rho),
                     capacity,
                     shortageCost
                     );
            }
            GSONUtility.<SKPMultinormal[]>saveInstanceToJSON(batch, batchFileName);
            break;
         }
      }  
   }
   
   public static SKPMultinormal[] retrieveBatch(String fileName) {
      SKPMultinormal[] instances = GSONUtility.<SKPMultinormal[]>retrieveJSONInstance(fileName, SKPMultinormal[].class);
      return instances;
   }
   
   /*
    * MILP
    */   
   
   public static void solveMILP(String fileName, int partitions, int simulationRuns, String folder) throws IloException {
      SKPMultinormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = folder+"/solved_multinormal_instances_MILP.json";
      SKPMultinormalMILPSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, partitions, simulationRuns);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPMultinormalMILPSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = folder+"/solved_multinormal_instances_MILP.csv";
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
            "instanceID, expectedValues, expectedWeights, covarianceWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, milpOptimalityGap, piecewisePartitions, "
            + "piecewiseSamples, milpMaxLinearizationError, simulatedLinearizationError,"
            + "cplexSolutionTimeMs, simplexIterations, exploredNodes\n";
      String body = "";
      
      for(SKPMultinormalMILPSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
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
         e.printStackTrace();
      }
   }
   
   private static SKPMultinormalMILPSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPMultinormalMILPSolvedInstance[] solvedInstances = GSONUtility.<SKPMultinormalMILPSolvedInstance[]>retrieveJSONInstance(fileName, SKPMultinormalMILPSolvedInstance[].class);
      return solvedInstances;
   }
   
   /*
    * DSKP
    */
   
   public static void solveDSKP(String fileName, String folder) {
      SKPMultinormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = folder+"/solved_multinormal_instances_DSKP.json";
      DSKPMultinormalSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved);
      System.out.println(GSONUtility.<DSKPMultinormalSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = folder+"/solved_multinormal_instances_DSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static DSKPMultinormalSolvedInstance[] solveBatchDSKP(SKPMultinormal[] instances, String fileName) {
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPMultinormalSolvedInstance>solved = new ArrayList<DSKPMultinormalSolvedInstance>();
      for(SKPMultinormal instance : instances) {
         solved.add(new DSKPMultinormal(instance, truncationQuantile).solve());
         System.out.println("Solved DSKP instance number "+solved.size());
         GSONUtility.<DSKPMultinormalSolvedInstance[]>saveInstanceToJSON(solved.toArray(new DSKPMultinormalSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPMultinormalSolvedInstance[solved.size()]);
   }

   private static void storeSolvedBatchToCSV(DSKPMultinormalSolvedInstance[] instances, String fileName) {
      String header = "instanceID, expectedValues, expectedWeights, covarianceWeights, capacity, shortageCost, solutionValue, solutionTimeMs, statesExplored\n";
      String body = "";
      
      for(DSKPMultinormalSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
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
         e.printStackTrace();
      }
   }
   
   private static DSKPMultinormalSolvedInstance[] retrieveSolvedBatchDSKP(String fileName) {
      DSKPMultinormalSolvedInstance[] solvedInstances = GSONUtility.<DSKPMultinormalSolvedInstance[]>retrieveJSONInstance(fileName, DSKPMultinormalSolvedInstance[].class);
      return solvedInstances;
   }
   
   protected static void storeBatchAsOPLDataFiles(SKPMultinormal[] instances, String OPLDataFileZipArchive, int partitions) {
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
