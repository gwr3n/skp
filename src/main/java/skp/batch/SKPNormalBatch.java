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
import skp.SKPNormal;
import skp.milp.SKPNormalMILP;
import skp.milp.SKPNormalMILPSolvedInstance;
import skp.sdp.DSKPNormal;
import skp.sdp.DSKPNormalSolvedInstance;
import skp.utililities.gson.GSONUtility;
import umontreal.ssj.probdist.UniformDist;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.randvar.UniformGen;
import umontreal.ssj.randvar.UniformIntGen;
import umontreal.ssj.rng.MRG32k3aL;

public class SKPNormalBatch {
   static final long[] seed = {1,2,3,4,5,6};
   static final private MRG32k3aL randGenerator = new MRG32k3aL();

   public static void main(String args[]) {
      String batchFileName = "scrap/normal_instances.json";
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
      SKPNormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "scrap/solvedNormalInstancesMILP.json";
      SKPNormalMILPSolvedInstance[] solvedBatch = solveBatchMILP(batch, fileNameSolved);
      
      solvedBatch = retrieveSolvedBatchMILP(fileNameSolved);
      System.out.println(GSONUtility.<SKPNormalMILPSolvedInstance[]>printInstanceAsGSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solvedInstancesMILP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
   }
   
   private static void storeSolvedBatchToCSV(SKPNormalMILPSolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValuesPerUnit, expectedWeights, stdWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "simulationRuns, milpSolutionValue, simulationRuns, milpOptimalityGap, piecewisePartitions, "
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
   
   private static SKPNormalMILPSolvedInstance[] solveBatchMILP(SKPNormal[] instances, String fileName) throws IloException {
      int partitions = 10;
      int simulationRuns = 100000;
      
      ArrayList<SKPNormalMILPSolvedInstance>solved = new ArrayList<SKPNormalMILPSolvedInstance>();
      for(SKPNormal instance : instances) {
         solved.add(new SKPNormalMILP(instance, partitions).solve(simulationRuns));
         GSONUtility.<SKPNormalMILPSolvedInstance[]>saveInstanceToGSON(solved.toArray(new SKPNormalMILPSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new SKPNormalMILPSolvedInstance[solved.size()]);
   }
   
   private static SKPNormalMILPSolvedInstance[] retrieveSolvedBatchMILP(String fileName) {
      SKPNormalMILPSolvedInstance[] solvedInstances = GSONUtility.<SKPNormalMILPSolvedInstance[]>retrieveInstance(fileName, SKPNormalMILPSolvedInstance[].class);
      return solvedInstances;
   }
   
   /*
    * DSKP
    */
   
   public static void solveDSKP(String fileName) {
      SKPNormal[] batch = retrieveBatch(fileName);
      
      String fileNameSolved = "scrap/solvedNormalInstancesDSKP.json";
      DSKPNormalSolvedInstance[] solvedBatch = solveBatchDSKP(batch, fileNameSolved);
      
      solvedBatch = retrieveSolvedBatchDSKP(fileNameSolved);
      System.out.println(GSONUtility.<DSKPNormalSolvedInstance[]>printInstanceAsGSON(solvedBatch));
      
      String fileNameSolvedCSV = "scrap/solvedInstancesDSKP.csv";
      storeSolvedBatchToCSV(solvedBatch, fileNameSolvedCSV);
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
   
   private static DSKPNormalSolvedInstance[] solveBatchDSKP(SKPNormal[] instances, String fileName) {
      double truncationQuantile = 0.999999999999999;
      ArrayList<DSKPNormalSolvedInstance>solved = new ArrayList<DSKPNormalSolvedInstance>();
      for(SKPNormal instance : instances) {
         solved.add(new DSKPNormal(instance, truncationQuantile).solve());
         GSONUtility.<DSKPNormalSolvedInstance[]>saveInstanceToGSON(solved.toArray(new DSKPNormalSolvedInstance[solved.size()]), fileName);
      }
      return solved.toArray(new DSKPNormalSolvedInstance[solved.size()]);
   }
   
   private static DSKPNormalSolvedInstance[] retrieveSolvedBatchDSKP(String fileName) {
      DSKPNormalSolvedInstance[] solvedInstances = GSONUtility.<DSKPNormalSolvedInstance[]>retrieveInstance(fileName, DSKPNormalSolvedInstance[].class);
      return solvedInstances;
   }   
   
   /**
    * Generate a batch of instances
    */
   
   public static void generateBatch(int numberOfInstances, int instanceSize, String fileName) {
      SKPNormal[] instances = SKPNormalBatch.generateInstances(numberOfInstances, instanceSize);
      GSONUtility.<SKPNormal[]>saveInstanceToGSON(instances, fileName);
   }
   
   public static SKPNormal[] retrieveBatch(String fileName) {
      SKPNormal[] instances = GSONUtility.<SKPNormal[]>retrieveInstance(fileName, SKPNormal[].class);
      return instances;
   }
   
   private static SKPNormal[] generateInstances(int numberOfInstances, int instanceSize){
      randGenerator.setSeed(seed);
      randGenerator.resetStartStream();
      SKPNormal[] instances = IntStream.iterate(0, i -> i + 1)
                                        .limit(numberOfInstances)
                                        .mapToObj(i -> new SKPNormal(
                                              (new RandomVariateGen(randGenerator, new UniformDist(0.1,10))).nextArrayOfDouble(instanceSize),
                                              (new RandomVariateGen(randGenerator, new UniformDist(15,70))).nextArrayOfDouble(instanceSize),
                                              UniformGen.nextDouble(randGenerator, 0.1, 0.5),
                                              UniformIntGen.nextInt(randGenerator, 100, 200),
                                              UniformGen.nextDouble(randGenerator, 50, 150)))
                                        .toArray(SKPNormal[]::new);
      return instances;
   }
}
