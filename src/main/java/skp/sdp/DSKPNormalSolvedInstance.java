package skp.sdp;

import skp.instance.SKPNormal;

public class DSKPNormalSolvedInstance {
   public SKPNormal instance;
   public double solutionValue;
   public long solutionTimeMs;
   public long statesExplored;
   
   public DSKPNormalSolvedInstance(SKPNormal instance, 
                                   double solutionValue,
                                   long solutionTimeMs,
                                   long statesExplored) {
      this.instance = instance;
      this.solutionValue = solutionValue;
      this.solutionTimeMs = solutionTimeMs;
      this.statesExplored = statesExplored;
   }
}
