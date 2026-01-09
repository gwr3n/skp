package skp.milp.instance;

import skp.instance.SKPNormal;

public class SKPNormalMILPSolvedInstance {
   public SKPNormal instance;
   public int[] optimalKnapsack;
   public double simulatedSolutionValue;
   public int simulationRuns;
   public double milpSolutionValue;
   public double milpOptimalityGap;
   public int piecewisePartitions;
   public int piecewiseSamples;
   public double milpMaxLinearizationError;
   public double simulatedLinearizationError;
   public double cplexSolutionTimeMs;
   public int simplexIterations;
   public int exploredNodes;
   
   public SKPNormalMILPSolvedInstance(
         SKPNormal instance,
         int[] optimalKnapsack,
         double simulatedSolutionValue,
         int simulationRuns,
         double milpSolutionValue,
         double milpOptimalityGap,
         int piecewisePartitions,
         int piecewiseSamples,
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
      this.piecewisePartitions = piecewisePartitions;
      this.piecewiseSamples = piecewiseSamples;
      this.milpMaxLinearizationError = milpMaxLinearizationError;
      this.simulatedLinearizationError = simulatedLinearizationError;
      this.cplexSolutionTimeMs = cplexSolutionTimeMs;
      this.simplexIterations = simplexIterations;
      this.exploredNodes = exploredNodes;
   }
   
}
