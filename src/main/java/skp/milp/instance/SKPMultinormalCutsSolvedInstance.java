package skp.milp.instance;

import skp.instance.SKPMultinormal;

public class SKPMultinormalCutsSolvedInstance {
   public SKPMultinormal instance;
   public int[] optimalKnapsack;
   public double simulatedSolutionValue;
   public int simulationRuns;
   public double milpSolutionValue;
   public double milpOptimalityGap;
   public int cuts;
   public double milpMaxLinearizationError;
   public double simulatedLinearizationError;
   public double cplexSolutionTimeMs;
   public int simplexIterations;
   public int exploredNodes;
   
   public SKPMultinormalCutsSolvedInstance(
         SKPMultinormal instance,
         int[] optimalKnapsack,
         double simulatedSolutionValue,
         int simulationRuns,
         double milpSolutionValue,
         double milpOptimalityGap,
         int cuts,
         double milpMaxLinearizationError,
         double simulatedLinearizationError,
         double cplexSolutionTimeMs,
         int simplexIterations,
         int exploredNodes) {
      this.instance = instance;
      this.optimalKnapsack = optimalKnapsack;
      this.simulatedSolutionValue = simulatedSolutionValue;
      this.simulationRuns = simulationRuns;
      this.milpSolutionValue = milpSolutionValue;
      this.milpOptimalityGap = milpOptimalityGap;
      this.cuts = cuts;
      this.milpMaxLinearizationError = milpMaxLinearizationError;
      this.simulatedLinearizationError = simulatedLinearizationError;
      this.cplexSolutionTimeMs = cplexSolutionTimeMs;
      this.simplexIterations = simplexIterations;
      this.exploredNodes = exploredNodes;
   }
}
