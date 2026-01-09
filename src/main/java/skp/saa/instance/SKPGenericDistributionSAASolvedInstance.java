package skp.saa.instance;

import skp.instance.SKPGenericDistribution;

public class SKPGenericDistributionSAASolvedInstance {
   public SKPGenericDistribution instance;
   public int[] optimalKnapsack;
   public double simulatedSolutionValue;
   public double solutionTimeMs;
   public double optGap1;
   public double optGap2;
   public int Nsmall;
   public int Nlarge;
   public int M;

   public SKPGenericDistributionSAASolvedInstance(SKPGenericDistribution instance, int[] optimalKnapsack, double simulatedSolutionValue, double solutionTimeMs, double optGap1, double optGap2, int Nsmall, int Nlarge, int M) {
       this.instance = instance;
       this.optimalKnapsack = optimalKnapsack;
       this.simulatedSolutionValue = simulatedSolutionValue;
       this.solutionTimeMs = solutionTimeMs;
       this.optGap1 = optGap1;
       this.optGap2 = optGap2;
       this.Nsmall = Nsmall;
       this.Nlarge = Nlarge;
       this.M = M;
   }
}
