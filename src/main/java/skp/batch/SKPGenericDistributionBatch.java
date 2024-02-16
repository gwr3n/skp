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
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

import ilog.concert.IloException;

import skp.instance.SKPGenericDistribution;
import skp.milp.SKPGenericDistributionMILP;
import skp.milp.instance.SKPGenericDistributionMILPSolvedInstance;
import skp.sdp.DSKPGenericDistribution;
import skp.sdp.instance.DSKPGenericDistributionSolvedInstance;
import skp.utililities.gson.GSONUtility;

import umontreal.ssj.probdist.DiscreteDistributionInt;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.probdist.UniformIntDist;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.randvar.RandomVariateGenInt;

public class SKPGenericDistributionBatch extends SKPBatch {
   
   public static void main(String args[]) {
      File folder = new File("scrap");
      if (!folder.exists()) {
        folder.mkdir();
      } 
      
      String batchFileName = "scrap/generic_distribution_instances.json";
      
      SKPGenericDistribution[] instances = generateInstances(batchFileName);
      
      int linearizationSamples = 100000;
      int simulationRuns = 100000;      
      try {
         solveMILP(instances, linearizationSamples, simulationRuns);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      solveDSKP(instances);
   }

   private static SKPGenericDistribution[] generateInstances(String batchFileName) {
      int instances = 10;
      int instanceSize = 10;
      
      Distribution expectedValuePerUnit = new UniformDist(0.1,10);
      Distribution expectedWeight = new UniformDist(15,70);
      Distribution coefficientOfVariation = new UniformDist(0.1, 0.5);
      DiscreteDistributionInt capacity = new UniformIntDist(100,200);
      Distribution shortageCost = new UniformDist(50,150);
      
      SKPGenericDistributionBatch batch = new SKPGenericDistributionBatch(expectedValuePerUnit, expectedWeight, coefficientOfVariation, capacity, shortageCost);
      return batch.generateBatch(instances, instanceSize, batchFileName);
   }
   
   protected Distribution coefficientOfVariation;
   
   public SKPGenericDistributionBatch(
         Distribution expectedValuePerUnit,
         Distribution expectedWeight,
         Distribution coefficientOfVariation,
         DiscreteDistributionInt capacity,
         Distribution shortageCost) {
      super(expectedValuePerUnit, expectedWeight, capacity, shortageCost);
      this.coefficientOfVariation = coefficientOfVariation;
   }
   
   /**
    * Generate a batch of instances
    */
   
   public SKPGenericDistribution[] generateBatch(int numberOfInstances, int instanceSize, String fileName) {
      SKPGenericDistribution[] instances = this.generateInstances(numberOfInstances, instanceSize);
      GSONUtility.<SKPGenericDistribution[]>saveInstanceToJSON(instances, fileName);
      return instances;
   }
   
   private SKPGenericDistribution[] generateInstances(int numberOfInstances, int instanceSize){
      randGenerator.setSeed(seed);
      randGenerator.resetStartStream();
      SKPGenericDistribution[] instances = IntStream.iterate(0, i -> i + 1)
                                                    .limit(numberOfInstances)
                                                    .mapToObj(i -> new SKPGenericDistribution(
                                                          (new RandomVariateGen(randGenerator, this.expectedValuePerUnit)).nextArrayOfDouble(instanceSize),
                                                          Arrays.stream((new RandomVariateGen(randGenerator, this.expectedWeight)).nextArrayOfDouble(instanceSize)).mapToObj(w -> new GammaDist(w, Math.sqrt(w/Math.pow((new RandomVariateGen(randGenerator, this.coefficientOfVariation)).nextDouble()*w,2)))).toArray(GammaDist[]::new),
                                                          (new RandomVariateGenInt(randGenerator, this.capacity)).nextInt(),
                                                          (new RandomVariateGen(randGenerator, this.shortageCost)).nextDouble()))
                                                    .toArray(SKPGenericDistribution[]::new);
      return instances;
   }
   
   public static void solveMILP(SKPGenericDistribution[] batch, int linearizationSamples, int simulationRuns) throws IloException {
      // SKPNormal[] batch = retrieveBatch(fileName); // Batch cannot be retrieved because Distribution[] is not Serializable
      
      String fileNameSolved = "scrap/solvedGenericDistributionInstancesMILP.json";
      SKPGenericDistributionMILPSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved, linearizationSamples, simulationRuns);
      
      // solvedBatch = retrieveSolvedBatchMILP(fileNameSolved); // Batch cannot be retrieved because Distribution[] is not Serializable
      System.out.println(GSONUtility.<SKPGenericDistributionMILPSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solvedGenericDistributionInstancesMILP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   public static SKPGenericDistribution[] retrieveBatch(String fileName) {
      SKPGenericDistribution[] instances = GSONUtility.<SKPGenericDistribution[]>retrieveJSONInstance(fileName, SKPGenericDistribution[].class);
      return instances;
   }
   
   private static SKPGenericDistributionMILPSolvedInstance[] solveBatchMILP(SKPGenericDistribution[] instances, String fileName, int linearizationSamples, int simulationRuns) throws IloException {
      ArrayList<SKPGenericDistributionMILPSolvedInstance>solved = new ArrayList<SKPGenericDistributionMILPSolvedInstance>();
      for(SKPGenericDistribution instance : instances) {
         solved.add(new SKPGenericDistributionMILP(instance, linearizationSamples).solve(simulationRuns));
         GSONUtility.<SKPGenericDistributionMILPSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPGenericDistributionMILPSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPGenericDistributionMILPSolvedInstance[solved.size()]);
   }
   
   private static SKPGenericDistributionMILPSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPGenericDistributionMILPSolvedInstance[] solvedInstances = GSONUtility.<SKPGenericDistributionMILPSolvedInstance[]>retrieveJSONInstance(fileName, SKPGenericDistributionMILPSolvedInstance[].class);
      return solvedInstances;
   }
   
   private static void storeSolvedBatchToCSV(SKPGenericDistributionMILPSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValuesPerUnit, expectedWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, milpOptimalityGap, cuts, "
            + "piecewiseSamples, milpMaxLinearizationError, simulatedLinearizationError,"
            + "cplexSolutionTimeMs, simplexIterations, exploredNodes\n";
      String body = "";
      
      for(SKPGenericDistributionMILPSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).map(d -> d.toString()).toArray()).replace(",", "\t")+ ", " +
                 s.instance.getCapacity()+ ", " +
                 s.instance.getShortageCost()+ ", " +
                 Arrays.toString(s.optimalKnapsack).replace(",", "\t")+ ", " +
                 s.simulatedSolutionValue + ", " +
                 s.simulationRuns + ", " +
                 s.milpSolutionValue + ", " +
                 s.milpOptimalityGap + ", " +
                 s.cuts + ", " +
                 s.piecewiseSamples  + ", " +
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
   
   /*
    * DSKP
    */
   
   public static void solveDSKP(SKPGenericDistribution[] batch) {
      // SKPNormal[] batch = retrieveBatch(fileName); // Batch cannot be retrieved because Distribution[] is not Serializable
      
      String fileNameSolved = "scrap/solvedGenericDistributionInstancesDSKP.json";
      DSKPGenericDistributionSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      // solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved); // Batch cannot be retrieved because Distribution[] is not Serializable
      System.out.println(GSONUtility.<DSKPGenericDistributionSolvedInstance[]>printInstanceAsJSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solvedGenericDistributionInstancesDSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static DSKPGenericDistributionSolvedInstance[] solveBatchDSKP(SKPGenericDistribution[] instances, String fileName) {
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPGenericDistributionSolvedInstance>solved = new ArrayList<DSKPGenericDistributionSolvedInstance>();
      for(SKPGenericDistribution instance : instances) {
         solved.add(new DSKPGenericDistribution(instance, truncationQuantile).solve());
         GSONUtility.<DSKPGenericDistributionSolvedInstance[]>saveInstanceToJSON(solved.toArray(new DSKPGenericDistributionSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPGenericDistributionSolvedInstance[solved.size()]);
   }

   private static void storeSolvedBatchToCSV(DSKPGenericDistributionSolvedInstance[] instances, String fileName) {
      String header = "instanceID, expectedValuesPerUnit, expectedWeights, capacity, shortageCost, solutionValue, solutionTimeMs, statesExplored\n";
      String body = "";
      
      for(DSKPGenericDistributionSolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValuesPerUnit()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).map(d -> d.toString()).toArray()).replace(",", "\t")+ ", " +
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
   
   private static DSKPGenericDistributionSolvedInstance[] retrieveSolvedBatchDSKP(String fileName) {
      DSKPGenericDistributionSolvedInstance[] solvedInstances = GSONUtility.<DSKPGenericDistributionSolvedInstance[]>retrieveJSONInstance(fileName, DSKPGenericDistributionSolvedInstance[].class);
      return solvedInstances;
   }
   
}
