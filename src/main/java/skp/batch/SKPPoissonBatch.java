/**
 * To run from Mac OS
 * 
 * -Djava.library.path=/Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
 * 
 * Environment variable
 * 
 * DYLD_LIBRARY_PATH /Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
 * 
 * @author gwren
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

import skp.instance.SKPPoisson;
import skp.milp.SKPPoissonMILP;
import skp.milp.SKPPoissonMILPSolvedInstance;
import skp.sdp.DSKPPoisson;
import skp.sdp.DSKPPoissonSolvedInstance;
import skp.utililities.gson.GSONUtility;

import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.randvar.UniformIntGen;
import umontreal.ssj.rng.MRG32k3aL;

public class SKPPoissonBatch {
   static final long[] seed = {1,2,3,4,5,6};
   static final private MRG32k3aL randGenerator = new MRG32k3aL();
   
   public static void main(String args[]) {
      String batchFileName = "scrap/poisson_instances.json";
      generateBatch(1, 100, batchFileName);
      
      try {
         solveMILP(batchFileName);
      } catch (IloException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      solveDSKP(batchFileName);
   }
   
   /*
    * MILP
    */
   
   public static void solveMILP(String fileName) throws IloException {
      SKPPoisson[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "scrap/solvedPoissonInstancesMILP.json";
      SKPPoissonMILPSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPPoissonMILPSolvedInstance[]>printInstanceAsGSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solvedPoissonInstancesMILP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static void storeSolvedBatchToCSV(SKPPoissonMILPSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValuesPerUnit, expectedWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, simulationRuns, milpOptimalityGap, piecewisePartitions, "
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
                 s.simulationRuns + ", " +
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
   
   private static SKPPoissonMILPSolvedInstance[] solveBatchMILP(SKPPoisson[] instances, String fileName) throws IloException {
      int partitions = 10;
      int linearizationSamples = 50000;
      int simulationRuns = 100000;
      
      ArrayList<SKPPoissonMILPSolvedInstance>solved = new ArrayList<SKPPoissonMILPSolvedInstance>();
      for(SKPPoisson instance : instances) {
         solved.add(new SKPPoissonMILP(instance, partitions, linearizationSamples).solve(simulationRuns));
         GSONUtility.<SKPPoissonMILPSolvedInstance[]>saveInstanceToGSON(solved.toArray(new SKPPoissonMILPSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPPoissonMILPSolvedInstance[solved.size()]);
   }
   
   private static SKPPoissonMILPSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPPoissonMILPSolvedInstance[] solvedInstances = GSONUtility.<SKPPoissonMILPSolvedInstance[]>retrieveInstance(fileName, SKPPoissonMILPSolvedInstance[].class);
      return solvedInstances;
   }
   
   /*
    * DSKP
    */
   
   public static void solveDSKP(String fileName) {
      SKPPoisson[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "scrap/solvedPoissonInstancesDSKP.json";
      DSKPPoissonSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved);
      System.out.println(GSONUtility.<DSKPPoissonSolvedInstance[]>printInstanceAsGSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solvedPoissonInstancesDSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
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
   
   private static DSKPPoissonSolvedInstance[] solveBatchDSKP(SKPPoisson[] instances, String fileName) {
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPPoissonSolvedInstance>solved = new ArrayList<DSKPPoissonSolvedInstance>();
      for(SKPPoisson instance : instances) {
         solved.add(new DSKPPoisson(instance, truncationQuantile).solve());
         GSONUtility.<DSKPPoissonSolvedInstance[]>saveInstanceToGSON(solved.toArray(new DSKPPoissonSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPPoissonSolvedInstance[solved.size()]);
   }
   
   private static DSKPPoissonSolvedInstance[] retrieveSolvedBatchDSKP(String fileName) {
      DSKPPoissonSolvedInstance[] solvedInstances = GSONUtility.<DSKPPoissonSolvedInstance[]>retrieveInstance(fileName, DSKPPoissonSolvedInstance[].class);
      return solvedInstances;
   }
   
   /**
    * Generate a batch of instances
    */
   
   public static void generateBatch(int numberOfInstances, int instanceSize, String fileName) {
      SKPPoisson[] instances = SKPPoissonBatch.generateInstances(numberOfInstances, instanceSize);
      GSONUtility.<SKPPoisson[]>saveInstanceToGSON(instances, fileName);
   }
   
   public static SKPPoisson[] retrieveBatch(String fileName) {
      SKPPoisson[] instances = GSONUtility.<SKPPoisson[]>retrieveInstance(fileName, SKPPoisson[].class);
      return instances;
   }
   
   private static SKPPoisson[] generateInstances(int numberOfInstances, int instanceSize){
      randGenerator.setSeed(seed);
      randGenerator.resetStartStream();
      SKPPoisson[] instances = IntStream.iterate(0, i -> i + 1)
                                        .limit(numberOfInstances)
                                        .mapToObj(i -> new SKPPoisson(
                                              (new RandomVariateGen(randGenerator, new UniformDist(0.1,10))).nextArrayOfDouble(instanceSize),
                                              (new RandomVariateGen(randGenerator, new UniformDist(15,70))).nextArrayOfDouble(instanceSize),
                                              UniformIntGen.nextInt(randGenerator, 100, 200),
                                              UniformGen.nextDouble(randGenerator, 50, 150)))
                                        .toArray(SKPPoisson[]::new);
      return instances;
   }   
}
