package skp.milp.instance;

import skp.instance.SKPGenericDistribution;

public class SKPScenarioBasedSolvedInstance {
   public SKPGenericDistribution instance;
   public int[] optimalKnapsack;
   public double simulatedSolutionValueMean;
   public double simulatedSolutionValueVariance;
   public int simulationRuns;
   public double milpSolutionValue;
   public double milpOptimalityGap;
   public double cplexSolutionTimeMs;
   public int simplexIterations;
   public int exploredNodes;
   
   public SKPScenarioBasedSolvedInstance(SKPGenericDistribution instance, int[] optimalKnapsack, double simulatedSolutionValueMean, double simulatedSolutionValueVariance, int simulationRuns, double milpSolutionValue, double milpOptimalityGap, double cplexSolutionTimeMs, int simplexIterations, int exploredNodes) {
      this.instance = instance;
      this.optimalKnapsack = optimalKnapsack;
      this.simulatedSolutionValueMean = simulatedSolutionValueMean;
      this.simulatedSolutionValueVariance = simulatedSolutionValueVariance;
      this.simulationRuns = simulationRuns;
      this.milpSolutionValue = milpSolutionValue;
      this.milpOptimalityGap = milpOptimalityGap;
      this.cplexSolutionTimeMs = cplexSolutionTimeMs;
      this.simplexIterations = simplexIterations;
      this.exploredNodes = exploredNodes;
   }
}
