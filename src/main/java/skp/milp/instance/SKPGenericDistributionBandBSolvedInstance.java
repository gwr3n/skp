package skp.milp.instance;

import skp.instance.SKPGenericDistribution;

public class SKPGenericDistributionBandBSolvedInstance {
   public SKPGenericDistribution instance;
   public int[] optimalKnapsack;
   public double simulatedSolutionValue;
   public int simulationRuns;
   public long solutionTimeMs;
   public int exploredNodes;
   public double optGap;
   
   public SKPGenericDistributionBandBSolvedInstance(
         SKPGenericDistribution instance,
         int[] optimalKnapsack,
         double simulatedSolutionValue,
         int simulationRuns,
         long solutionTimeMs,
         int exploredNodes,
         double optGap) {
      this.instance = instance;
      this.optimalKnapsack = optimalKnapsack;
      this.simulatedSolutionValue = simulatedSolutionValue;
      this.simulationRuns = simulationRuns;
      this.solutionTimeMs = solutionTimeMs;
      this.exploredNodes = exploredNodes;
      this.optGap = optGap;
   }
}
