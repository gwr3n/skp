package skp.sdp.instance;

import skp.instance.SKPMultinormal;

public class DSKPMultinormalSolvedInstance {
   public SKPMultinormal instance;
   public double solutionValue;
   public long solutionTimeMs;
   public long statesExplored;
   
   public DSKPMultinormalSolvedInstance(SKPMultinormal instance,
                                        double solutionValue,
                                        long solutionTimeMs,
                                        long statesExplored) {
      this.instance = instance;
      this.solutionValue = solutionValue;
      this.solutionTimeMs = solutionTimeMs;
      this.statesExplored = statesExplored;
   }

}
