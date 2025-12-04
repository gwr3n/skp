package skp.saa.instance;

import skp.instance.SKPGenericDistribution;

public class SKPGenericDistributionSAA_LDSolvedInstance {
   public SKPGenericDistribution instance;
   
   public int[] optimalKnapsack;
   public double simulatedSolutionValue;
   
   // Phase 0 (LD sizing)
   public long N_LD_initial, N_LD_final, N_LD_max;
   public long N_start, N_last, N_attempts;
   public String phase0StopReason;
   
   // Phase 1 (replication loop at the chosen N)
   public long N_phase1, Mfinal;
   public double vBar, gHat, relGap1, relGap2;
   public String phase1StopReason;
   public long solutionTimeMs;
   
   public SKPGenericDistributionSAA_LDSolvedInstance(
         SKPGenericDistribution instance, 
         int[] optimalKnapsack, 
         double simulatedSolutionValue,
         long N_LD_initial,
         long N_LD_final,
         long N_LD_max,
         long N_start,
         long N_last,
         long N_attempts,
         String phase0StopReason,
         long N_phase1,
         long Mfinal,
         double vBar,
         double gHat,
         double relGap1,
         double relGap2,
         String phase1StopReason,
         long solutionTimeMs
         ) {
       this.instance = instance;
       this.optimalKnapsack = optimalKnapsack;
       this.simulatedSolutionValue = simulatedSolutionValue;
       this.solutionTimeMs = solutionTimeMs;
       this.N_LD_initial = N_LD_initial;
       this.N_LD_final = N_LD_final;
       this.N_LD_max = N_LD_max;
       this.N_start = N_start;
       this.N_last = N_last;
       this.N_attempts = N_attempts;
       this.phase0StopReason = phase0StopReason;
       this.N_phase1 = N_phase1;
       this.Mfinal = Mfinal;
       this.vBar = vBar;
       this.gHat = gHat;
       this.relGap1 = relGap1;
       this.relGap2 = relGap2;
       this.phase1StopReason = phase1StopReason;
       this.solutionTimeMs = solutionTimeMs;
   }
}
