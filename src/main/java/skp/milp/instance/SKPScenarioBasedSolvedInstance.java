package skp.milp.instance;

import skp.instance.SKPGenericDistribution;

public class SKPScenarioBasedSolvedInstance {
   public SKPGenericDistribution instance;
   public int[] optimalKnapsack;
   public double simulatedSolutionValue;
   public int simulationRuns;
   public double milpSolutionValue;
   public double milpOptimalityGap;
   public double cplexSolutionTimeMs;
   public int simplexIterations;
   public int exploredNodes;
   
   public SKPScenarioBasedSolvedInstance(SKPGenericDistribution instance, int[] optimalKnapsack, double simulatedSolutionValue, int simulationRuns, double milpSolutionValue, double milpOptimalityGap, double cplexSolutionTimeMs, int simplexIterations, int exploredNodes) {
      this.instance = instance;
      this.optimalKnapsack = optimalKnapsack;
      this.simulatedSolutionValue = simulatedSolutionValue;
      this.simulationRuns = simulationRuns;
      this.milpSolutionValue = milpSolutionValue;
      this.milpOptimalityGap = milpOptimalityGap;
      this.cplexSolutionTimeMs = cplexSolutionTimeMs;
      this.simplexIterations = simplexIterations;
      this.exploredNodes = exploredNodes;
   }
}
