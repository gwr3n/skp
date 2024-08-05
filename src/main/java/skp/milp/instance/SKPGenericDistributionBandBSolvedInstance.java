package skp.milp.instance;

import skp.instance.SKPGenericDistribution;

public class SKPGenericDistributionBandBSolvedInstance {
   public SKPGenericDistribution instance;
   public int[] optimalKnapsack;
   public double simulatedSolutionValue;
   public int simulationRuns;
   public long solutionTimeMs;
   public int exploredNodes;
   
   public SKPGenericDistributionBandBSolvedInstance(
         SKPGenericDistribution instance,
         int[] optimalKnapsack,
         double simulatedSolutionValue,
         int simulationRuns,
         long solutionTimeMs,
         int exploredNodes) {
      this.instance = instance;
      this.optimalKnapsack = optimalKnapsack;
      this.simulatedSolutionValue = simulatedSolutionValue;
      this.simulationRuns = simulationRuns;
      this.solutionTimeMs = solutionTimeMs;
      this.exploredNodes = exploredNodes;
   }
}
