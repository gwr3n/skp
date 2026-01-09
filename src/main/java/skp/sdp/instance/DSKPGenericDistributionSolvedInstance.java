package skp.sdp.instance;

import skp.instance.SKPGenericDistribution;

public class DSKPGenericDistributionSolvedInstance{
   public SKPGenericDistribution instance;
   public double solutionValue;
   public long solutionTimeMs;
   public long statesExplored;
   
   public DSKPGenericDistributionSolvedInstance(SKPGenericDistribution instance,
                                                double solutionValue,
                                                long solutionTimeMs,
                                                long statesExplored) {
      this.instance = instance;
      this.solutionValue = solutionValue;
      this.solutionTimeMs = solutionTimeMs;
      this.statesExplored = statesExplored;
   }
}
