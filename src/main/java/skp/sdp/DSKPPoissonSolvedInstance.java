package skp.sdp;

import skp.instance.SKPPoisson;

public class DSKPPoissonSolvedInstance {
   public SKPPoisson instance;
   public double solutionValue;
   public long solutionTimeMs;
   public long statesExplored;
   
   public DSKPPoissonSolvedInstance(SKPPoisson instance, 
                                   double solutionValue,
                                   long solutionTimeMs,
                                   long statesExplored) {
      this.instance = instance;
      this.solutionValue = solutionValue;
      this.solutionTimeMs = solutionTimeMs;
      this.statesExplored = statesExplored;
   }
}
