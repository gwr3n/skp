package skp.milp;

import skp.instance.SKPDistribution;
import skp.milp.instance.SKPDistributionMILPSolvedInstance;
import skp.utililities.gson.GSONUtility;

public class SKPDistributionMILP{
   int partitions;
   int linearizationSamples;
   
   int[] optimalKnapsack;
   double milpSolutionValue;
   double milpOptimalityGap;
   double milpMaxLinearizationError;
   
   SKPDistribution instance;
   
   public SKPDistributionMILP(SKPDistribution instance, int partitions, int linearizationSamples){
      this.instance = instance;
      this.partitions = partitions;
      this.linearizationSamples = linearizationSamples;
   }
   
   double gurobiSolutionTimeMs;
   int simplexIterations;
   int exploredNodes;

   public int[] getOptimalKnapsack() {
      return this.optimalKnapsack;
   }
   
   public double getMILPSolutionValue() {
      return this.milpSolutionValue;
   }
   
   public double getMILPOptimalityGap() {
      return this.milpOptimalityGap;
   }
   
   public double getMILPMaxLinearizationError() {
      return this.milpMaxLinearizationError;
   }

   void computeMILPMaxLinearizationError() {
      // TODO Auto-generated method stub
      // use PiecewiseFirstOrderLossFunction.getMaxApproximationError 
      // by feeding the optimal knapsack weight distributions to the constructor.
      System.err.print("Unimplemented method");
   }
   
   public SKPDistributionMILPSolvedInstance solve(int simulationRuns) {
      return null;
   }
   
   public static void main(String args[]) {

      SKPDistribution instance = SKPDistribution.getTestInstance();
      
      try {
         SKPDistributionMILP sskp = new SKPDistributionMILP(instance, 10, 100000);
         
         System.out.println(GSONUtility.<SKPDistributionMILPSolvedInstance>printInstanceAsJSON(sskp.solve(100000)));
      } catch (Exception e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }
}
