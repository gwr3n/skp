package skp.batch;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;

import ilog.concert.IloException;

import skp.instance.SKPGenericDistribution;
import skp.saa.SKPGenericDistributionSAA;
import skp.saa.instance.SKPGenericDistributionSAASolvedInstance;
import skp.utilities.gson.GSONUtility;

public class SKPGenericDistributionSAABatch {
   static SKPGenericDistributionSAASolvedInstance[] solveBatchMILPIterativeCuts(SKPGenericDistribution[] instances, String fileName, int Nsmall, int Nlarge, int M, double tolerance) throws IloException {
      /*
       * Sequential
       *
      ArrayList<SKPGenericDistributionSAASolvedInstance>solved = new ArrayList<SKPGenericDistributionSAASolvedInstance>();
      for(SKPGenericDistribution instance : instances) {
         solved.add(new SKPGenericDistributionSAA(instance).solve(Nsmall, Nlarge, M, tolerance));
         
      }
      GSONUtility.<SKPGenericDistributionCutsSolvedInstance[]>saveInstanceToJSON(solved.toArray(new SKPGenericDistributionCutsSolvedInstance[solved.size()]), fileName);
      return solved.toArray(new SKPGenericDistributionSAASolvedInstance[solved.size()]);
      */
      
      /*
       * Parallel
       */
      SKPGenericDistributionSAASolvedInstance[] solved = Arrays.stream(instances)
                                                               .parallel()
                                                               .map(instance -> {
           return new SKPGenericDistributionSAA(instance).solve(Nsmall, Nlarge, M, tolerance);
      }).toArray(SKPGenericDistributionSAASolvedInstance[]::new);
      GSONUtility.<SKPGenericDistributionSAASolvedInstance[]>saveInstanceToJSON(solved, fileName);
      return solved;
   }
   
   static void storeSolvedBatchToCSV(SKPGenericDistributionSAASolvedInstance[] instances, String fileName) {
      String header = 
            "instanceID, expectedValues, expectedWeights, "
            + "capacity, shortageCost, optimalKnapsack, simulatedSolutionValue, "
            + "solutionTimeMs, optGap1, optGap2, N, N', M\n";
      String body = "";
      
      for(SKPGenericDistributionSAASolvedInstance s : instances) {
         body += s.instance.getInstanceID() + ", " +
                 Arrays.toString(s.instance.getExpectedValues()).replace(",", "\t")+ ", " +
                 Arrays.toString(Arrays.stream(s.instance.getWeights()).map(d -> d.toString()).toArray()).replace(",", "\t")+ ", " +
                 s.instance.getCapacity()+ ", " +
                 s.instance.getShortageCost()+ ", " +
                 Arrays.toString(s.optimalKnapsack).replace(",", "\t")+ ", " +
                 s.simulatedSolutionValue + ", " +
                 s.solutionTimeMs + ", " +
                 s.optGap1 + ", " +
                 s.optGap2 + ", " +
                 s.Nsmall + ", " +
                 s.Nlarge + ", " +
                 s.M + "\n";
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
}
