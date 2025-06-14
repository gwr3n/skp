package skp.milp.instance;

import skp.instance.SKPMultinormal;

public class SKPMultinormalScenarioBasedSolvedInstance {
   public SKPMultinormal instance;
   public int[] optimalKnapsack;
   public double simulatedSolutionValueMean1;
   public double simulatedSolutionValueMean2;
   public double simulatedSolutionValueVariance1;
   public double simulatedSolutionValueVariance2;
   public int simulationRuns;
   public double milpSolutionValue;
   public double milpOptimalityGap;
   public double cplexSolutionTimeMs;
   public int simplexIterations;
   public int exploredNodes;
   
   public SKPMultinormalScenarioBasedSolvedInstance(SKPMultinormal instance, int[] optimalKnapsack, double simulatedSolutionValueMean1, double simulatedSolutionValueMean2, double simulatedSolutionValueVariance1, double simulatedSolutionValueVariance2, int simulationRuns, double milpSolutionValue, double milpOptimalityGap, double cplexSolutionTimeMs, int simplexIterations, int exploredNodes) {
      this.instance = instance;
      this.optimalKnapsack = optimalKnapsack;
      this.simulatedSolutionValueMean1 = simulatedSolutionValueMean1;
      this.simulatedSolutionValueMean2 = simulatedSolutionValueMean2;
      this.simulatedSolutionValueVariance1 = simulatedSolutionValueVariance1;
      this.simulatedSolutionValueVariance2 = simulatedSolutionValueVariance2;
      this.simulationRuns = simulationRuns;
      this.milpSolutionValue = milpSolutionValue;
      this.milpOptimalityGap = milpOptimalityGap;
      this.cplexSolutionTimeMs = cplexSolutionTimeMs;
      this.simplexIterations = simplexIterations;
      this.exploredNodes = exploredNodes;
   }
}
