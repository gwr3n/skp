package skp.saa;

import java.util.ArrayList;
import java.util.List;

import ilog.concert.IloException;
import skp.instance.SKPGenericDistribution;
import skp.milp.SKPGenericDistributionScenarioBased;
import skp.milp.instance.SKPGenericDistributionScenarioBasedSolvedInstance;
import skp.saa.instance.SKPGenericDistributionSAA_LDSolvedInstance;
import skp.sim.SimulateGenericDistribution;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

/**  SAA + LD rule (Section 2) + two CLT gap estimators (Section 3). */
public final class SKPGenericDistributionSAA_LD {

    /* ---------------- user parameters ---------------- */
    private static final double alpha  = 0.05;      // LD tail
    private static final double delta  = 0.0;       // usually 0
    private static final int    N0     = 64;        // start N
    private static final int    Mmax   = 1000;      // maximum replications
    private static final int    Nprime = 100_000;   // evaluation sample
    private static       int    NCap   = Nprime;    // Maximum N
    private static final int    warmUp = 32;        // for CLT
    private static final double relTol = 1e-4;      // relative gap
    private static final long   wallMs = 10*60_000L; // time limit in milliseconds
    private static final long   modMs  = wallMs/33;  // time limit in milliseconds to solve 
                                                     // each SAA MILP (we need to solve at 
                                                     // least 33 problems for CLT)

    /* ---------------- members ---------------- */
    private final SKPGenericDistribution inst;
    private final RandomStream           rng;

    private final List<SKPGenericDistributionScenarioBasedSolvedInstance> reps = new ArrayList<>();
    private final List<double[][]> smallScen = new ArrayList<>(); // size-N scenarios

    private double sigma2Max = 0.0;
    private double sumV      = 0.0;
    private int    bestIdx   = -1;
    //private double bestMean1 = Double.NEGATIVE_INFINITY;
    //private double bestAbs   = 0.0;     // |Mean1|
    private double bestMean2 = Double.NEGATIVE_INFINITY;
    private double bestAbs2  = 0.0;       // |Mean2|    

    public SKPGenericDistributionSAA_LD(SKPGenericDistribution inst) {
        this.inst = inst;
        MRG32k3aL r = new MRG32k3aL();
        r.setSeed(new long[]{12345,24513,24531,42531,35124,32451});
        this.rng = r;
    }

    /*==================================================================*/
    public SKPGenericDistributionSAA_LDSolvedInstance solve() {
       
        long wall0 = System.currentTimeMillis();
        int  N     = N0;

        /* ---------- Phase 0 : enlarge N until LD requirement ---------- */
        class Phase0MiniLog {
           long nldFirst=-1, nldFinal=-1, nldMax=-1; 
           long N_start=N0, N_last=N0, N_attempts=1;
           String stop; // "SUCCESS", "MMAX", "TIMEOUT", "NCAP"
        }
        Phase0MiniLog p0 = new Phase0MiniLog();
        while (true) {
            
            if (System.currentTimeMillis() - wall0 > wallMs) {
               p0.stop  = "TIMEOUT";
               p0.N_last = N;
               break;
            }

            try {
               makeReplication(N);
            } catch (IloException e) {
               // TODO Auto-generated catch block
               e.printStackTrace();
            }

            if (reps.size() >= Mmax) {
               p0.stop = "MMAX";
               p0.N_last = N;
               System.out.println("DEBUG LD: reached Mmax");
               break;
            }
            
            if (sigma2Max == 0.0) {
               // Need enough evidence before declaring degeneracy
               if (reps.size() >= warmUp && bestIdx >= 0) {
                   // 1) Check all knapsacks identical to incumbent
                   boolean allSame = true;
                   int[] incX = reps.get(bestIdx).optimalKnapsack;
                   for (int m = 0; m < reps.size(); m++) {
                       int[] x = reps.get(m).optimalKnapsack;
                       if (x.length != incX.length) { allSame = false; break; }
                       for (int i = 0; i < x.length; i++) {
                           if (x[i] != incX[i]) { allSame = false; break; }
                       }
                       if (!allSame) break;
                   }

                   // 2) Confirm all pairwise diff variances are zero against incumbent
                   boolean allZeroVar = true;
                   if (allSame) {
                       for (int m = 0; m < reps.size(); m++) {
                           double var = diffVariance(incX, reps.get(m).optimalKnapsack, smallScen.get(m));
                           if (var > 0.0) { allZeroVar = false; break; }
                       }
                   }

                   if (allSame && allZeroVar) {
                       // True degeneracy: N_LD = 0 for any N
                       if (p0.nldFirst < 0) p0.nldFirst = 0;
                       p0.nldMax   = Math.max(p0.nldMax, 0);
                       p0.nldFinal = 0;
                       p0.stop     = "SUCCESS_SIGMA2ZERO";
                       p0.N_last   = N;
                       System.out.println("DEBUG LD: σ̂²max=0 with stable incumbent and zero diff–variance; N_LD=0");
                       break;
                   }
               }
               // Otherwise keep sampling to collect more evidence
               continue;
           }

            double epsAbs = Math.max(relTol * bestAbs2, 1e-8);    // ε relative
            int    k      = inst.getItems();
            double logS   = k * Math.log(2.0);                   // |S| = 2^k
            double gamma  = (epsAbs-delta)*(epsAbs-delta)
                          / (3.0 * sigma2Max);
            int    Nld    = (int)Math.ceil( (logS - Math.log(alpha)) / gamma );

            System.out.printf("DEBUG LD:  N=%d σ̂²max=%.3g γ̂=%.3g N_LD=%d%n",
                              N, sigma2Max, gamma, Nld);
            
            /* ------- logging ----------------------- */
            if (sigma2Max > 0 && p0.nldFirst < 0) 
               p0.nldFirst = Nld; 
            p0.nldMax = Math.max(p0.nldMax, Nld);
            p0.nldFinal = Nld; // overwrite each time; it will hold the last value at exit.
            /* ------- end logging ------------------- */

            if (N >= Nld) {
               p0.stop = "SUCCESS";
               p0.N_last = N;
               System.out.println("DEBUG LD: satisfied N_LD requirement");
               break;
            }
            
            long nextN = ((long) N) << 1; // or: long nextN = 2L * N;

            if (nextN > NCap) {
                if (N < NCap) {
                    System.out.println("DEBUG LD: capping N at NCap=" + NCap);
                    p0.N_attempts++;
                    N = NCap;
                    p0.N_last = N;
                    reset();
                    continue; // try with NCap before deciding
                } else {
                    System.out.println("DEBUG LD: cannot increase N beyond NCap; LD may be unattainable");
                    p0.stop  = "NCAP";
                    p0.N_last = N;
                    break;
                }
            }
            
            System.out.println("DEBUG LD: increasing N to " + nextN);
            p0.N_attempts++;
            N <<= 1;
            p0.N_last = N;
            reset();
        }

        /* ---------- Phase 1 : CLT gap loop (both estimators) ---------- */
        class Phase1MiniLog {
           int    N, Mfinal;
           double vBar = Double.NaN, gHat = Double.NaN, relGap1 = Double.NaN, relGap2 = Double.NaN;   // gap2 can be NaN if unused
           long   wallMs;
           String stop;                     // "CERTIFIED", "MMAX", "TIMEOUT"
        }
        Phase1MiniLog p1 = new Phase1MiniLog();
        p1.N      = N;
        
        double finalGap1Rel=Double.POSITIVE_INFINITY, finalGap2Rel=Double.POSITIVE_INFINITY;
        double finalEvalMean = Double.NaN;

        while (true) {

            if (reps.size() < warmUp) { 
               try {
               makeReplication(N);
               } catch (IloException e) {
                  // TODO Auto-generated catch block
                  e.printStackTrace();
               } 
               continue; 
            }

            int    M = reps.size();
            double vBar = sumV / M;
            double s2_M = varianceV(vBar); // This is Var_hat(vBar)

            /*----- estimator 1 : fresh evaluation of incumbent ----------*/
            double[] eval = freshEvaluation(Nprime);
            double   gHat = eval[0];
            double   varN = eval[1] / Nprime;

            double point1 = vBar - gHat;
            double varTot = varN + s2_M; 
            double z      = NormalDist.inverseF01(1.0 - alpha); 
            double gap1   = point1 + z * Math.sqrt(varTot);

            /*----- estimator 2 : uses small-scenario samples ------------*/
            double[] gms = new double[M];
            double gBarN = 0.0;
            for (int m = 0; m < M; m++) {
                gms[m] = computeSampleAverage(reps.get(bestIdx).optimalKnapsack, smallScen.get(m));
                gBarN += gms[m];
            }
            gBarN /= M;

            double diffMean = vBar - gBarN;          // point estimator 2

            double s2bar = 0.0;
            for (int m = 0; m < M; m++) {
                double vm = reps.get(m).milpSolutionValue;
                s2bar += Math.pow((vm - gms[m]) - diffMean, 2);
            }
            s2bar = s2bar / (M * (M - 1));               // \bar S²_M / M
            double gap2 = diffMean + z * Math.sqrt(s2bar);
            
            /* ------- logging ----------------------- */
            double abs_gHat = Math.max(1e-12, gHat);
            p1.vBar = vBar;   // \bar v_N
            p1.gHat = gHat;   // \hat g_{N'}(\hat x)
            p1.relGap1 = gap1/Math.max(1e-12, gHat);
            p1.relGap2 = gap2/Math.max(1e-12, gHat);   
            /* ------- end logging ------------------- */
            
            System.out.printf(
              "DEBUG GAP: rep=%2d  ĝ=%.6f  gap1=%.3e  gap2=%.3e  rel1=%.3e%n rel2=%.3e%n",
              M, gHat, gap1, gap2, gap1/abs_gHat, gap2/abs_gHat);
            
            if (p1.relGap1 < relTol /* && p1.relGap2 < relTol */) {
                finalGap1Rel = p1.relGap1;
                finalGap2Rel = p1.relGap2;
                finalEvalMean= gHat;
                System.out.println("DEBUG: certified optimality");
                
                p1.stop = "CERTIFIED";
                break;
            }
            if (M >= Mmax) {
                finalGap1Rel = p1.relGap1;
                finalGap2Rel = p1.relGap2;
                finalEvalMean= gHat;
                System.out.println("DEBUG: reached Mmax");
                
                p1.stop = "MMAX";
                break;
            }
            if (System.currentTimeMillis()-wall0 > wallMs) {
               finalGap1Rel = p1.relGap1;
               finalGap2Rel = p1.relGap2;
               finalEvalMean= gHat;
               System.out.println("DEBUG: timeout");
                
               p1.stop = "TIMEOUT";
               break;
            }
            try {
               makeReplication(N);
            } catch (IloException e) {
               // TODO Auto-generated catch block
               e.printStackTrace();
            } 
        }
        long sec = (System.currentTimeMillis()-wall0)/1000;
        
        /* ------- logging ----------------------- */
        p1.Mfinal = reps.size();
        p1.wallMs = sec*1000L;
        /* ------- end logging ------------------- */

        System.out.printf("TERMINATE: relGap1 %.3e  relGap2 %.3e  "
                          +"rep=%d  N=%d  time=%ds%n",
                          finalGap1Rel, finalGap2Rel, reps.size(), N, sec);

        /* -------- build solved-instance --------------------------------*/
        SKPGenericDistributionScenarioBasedSolvedInstance inc = reps.get(bestIdx);
        return new SKPGenericDistributionSAA_LDSolvedInstance(
              inst,
              inc.optimalKnapsack,
              finalEvalMean,
              p0.nldFirst,
              p0.nldFinal,
              p0.nldMax,
              p0.N_start,
              p0.N_last,
              p0.N_attempts,
              p0.stop,
              p1.N,
              p1.Mfinal,
              p1.vBar,
              p1.gHat,
              p1.relGap1,
              p1.relGap2,
              p1.stop,
              p1.wallMs
              );
    }

    /*==================== helpers ====================================*/
    private void makeReplication(int N) throws IloException {

       SKPGenericDistributionScenarioBased saa =
              new SKPGenericDistributionScenarioBased(inst, N, rng);
       SKPGenericDistributionScenarioBasedSolvedInstance sol = saa.solve(Nprime);
       
       if(sol.cplexSolutionTimeMs > modMs) {
           /*System.out.printf("DEBUG: MILP solve time %fms exceeded limit %dms%n",
                             sol.cplexSolutionTimeMs, modMs);*/
           NCap = N; // force capping N next time
       }

       reps.add(sol);
       smallScen.add(saa.scenarios);           // keep scenarios for diff–variance
       sumV += sol.milpSolutionValue;

       /* --------- update σ̂²max with variance of difference ------------- */
       if (bestIdx >= 0) {                     // incumbent already defined
           double diffVar = diffVariance(
                   reps.get(bestIdx).optimalKnapsack,
                   sol.optimalKnapsack,
                   saa.scenarios);
           sigma2Max = Math.max(sigma2Max, diffVar);
       }

       // incumbent update on Mean2 (align with SAA)
       if (sol.simulatedSolutionValueMean2 > bestMean2) {
           bestMean2 = sol.simulatedSolutionValueMean2;
           bestAbs2  = Math.abs(bestMean2);
           bestIdx   = reps.size() - 1;
           System.out.printf("DEBUG INC: new incumbent at rep %d Mean2=%.6f%n",
                             reps.size(), bestMean2);

           /* incumbent changed ⇒ recompute σ̂²max against *all* replications */
           sigma2Max = 0.0;
           for (int m = 0; m < reps.size(); m++) {
               double var = diffVariance(
                       reps.get(bestIdx).optimalKnapsack,
                       reps.get(m).optimalKnapsack,
                       smallScen.get(m));
               sigma2Max = Math.max(sigma2Max, var);
           }
       }
   }


    private double[] freshEvaluation(int nSim) {
        SimulateGenericDistribution sim = new SimulateGenericDistribution(inst);
        return sim.simulateMeanVariance(
                  reps.get(bestIdx).optimalKnapsack, nSim, rng);
    }

    private double varianceV(double vBar) {
       double s = 0.0;
       for (int i = 0; i < reps.size(); i++) {
           SKPGenericDistributionScenarioBasedSolvedInstance r = reps.get(i);
           s += (r.milpSolutionValue - vBar) * (r.milpSolutionValue - vBar);
       }
       return s / (reps.size() * (reps.size() - 1));
   }

    private double computeSampleAverage(int[] knap, double[][] scen) {
       double det = 0.0;
       for (int i = 0; i < knap.length; i++)
           if (knap[i] == 1) det += inst.getExpectedValues()[i]; // deterministic part

       double sum = 0.0;
       for (double[] row : scen) {
           double wt = 0.0;
           for (int i = 0; i < knap.length; i++)
               if (knap[i] == 1) wt += row[i];
           double penalty = inst.getShortageCost() * Math.max(0, wt - inst.getCapacity());
           sum += (det - penalty);
       }
       return sum / scen.length;
   }

    private void reset() {
        reps.clear();
        smallScen.clear();
        sigma2Max=0; sumV=0; bestIdx=-1;
        bestMean2=Double.NEGATIVE_INFINITY; bestAbs2=0;
        rng.resetStartStream();
    }
    
    /* ---------- new helper: variance of incumbent – otherSolution over N scenarios */
    private double diffVariance(int[] inc, int[] x,
                                double[][] scen) {
        int N = scen.length;
        double[] diff = new double[N];

        for (int s = 0; s < N; s++) {
            double vInc = 0.0, vX = 0.0, wtInc = 0.0, wtX = 0.0;

            /* deterministic part */
            for (int i = 0; i < inc.length; i++) {
                if (inc[i]==1) vInc += inst.getExpectedValues()[i];
                if (x  [i]==1) vX   += inst.getExpectedValues()[i];
            }
            /* scenario–dependant weight */
            for (int i = 0; i < inc.length; i++) {
                if (inc[i]==1) wtInc += scen[s][i];
                if (x  [i]==1) wtX   += scen[s][i];
            }
            vInc -= inst.getShortageCost() * Math.max(0, wtInc - inst.getCapacity());
            vX   -= inst.getShortageCost() * Math.max(0, wtX   - inst.getCapacity());
            diff[s] = vInc - vX;
        }
        double mean = 0.0;
        for (double d : diff) mean += d;
        mean /= N;

        double var = 0.0;
        for (double d : diff) var += (d-mean)*(d-mean);
        return var / Math.max(1, N - 1);                     // unbiased N-denominator
    }


    /* ---------------- quick test ------------------------------------ */
    public static void main(String[] args) throws IloException {
        SKPGenericDistribution inst =
              SKPGenericDistribution.getTestInstance();

        SKPGenericDistributionSAA_LD solver =
              new SKPGenericDistributionSAA_LD(inst);

        SKPGenericDistributionSAA_LDSolvedInstance res = solver.solve();
        System.out.println("\nJSON RESULT:\n"+
                           GSONUtility.printInstanceAsJSON(res));
    }
}
