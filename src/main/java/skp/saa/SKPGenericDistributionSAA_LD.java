package skp.saa;

import java.util.ArrayList;
import java.util.List;

import ilog.concert.IloException;
import skp.instance.SKPGenericDistribution;
import skp.milp.SKPGenericDistributionScenarioBased;
import skp.milp.instance.SKPGenericDistributionScenarioBasedSolvedInstance;
import skp.saa.instance.SKPGenericDistributionSAASolvedInstance;
import skp.sim.SimulateGenericDistribution;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

/**  SAA + LD rule (Section 2) + two CLT gap estimators (Section 3). */
public final class SKPGenericDistributionSAA_LD {

    /* ---------------- user parameters ---------------- */
    private static final double alpha  = 0.05;   // LD tail
    private static final double delta  = 0.0;    // usually 0
    private static final int    N0     = 64;     // start N
    private static final int    Mmax   = 100;    // maximum replications
    private static final int    Nprime = 10_000; // evaluation sample
    private static final int    warmUp = 32;     // for CLT
    private static final double relTol = 1e-4;   // relative gap
    private static final long   wallMs = 10*60_000L;

    /* ---------------- members ---------------- */
    private final SKPGenericDistribution inst;
    private final RandomStream           rng;

    private final List<SKPGenericDistributionScenarioBasedSolvedInstance> reps
            = new ArrayList<>();
    private final List<double[][]> smallScen = new ArrayList<>(); // size-N scenarios

    private double sigma2Max = 0.0;
    private double sumV      = 0.0;
    private int    bestIdx   = -1;
    private double bestMean1 = Double.NEGATIVE_INFINITY;
    private double bestAbs   = 0.0;     // |Mean1|

    public SKPGenericDistributionSAA_LD(SKPGenericDistribution inst) {
        this.inst = inst;
        MRG32k3aL r = new MRG32k3aL();
        r.setSeed(new long[]{12345,24513,24531,42531,35124,32451});
        this.rng = r;
    }

    /*==================================================================*/
    public SKPGenericDistributionSAASolvedInstance solve() {

        long wall0 = System.currentTimeMillis();
        int  N     = N0;

        /* ---------- Phase 0 : enlarge N until LD requirement ---------- */
        while (true) {

            try {
               makeReplication(N);
            } catch (IloException e) {
               // TODO Auto-generated catch block
               e.printStackTrace();
            }

            if (reps.size() >= Mmax) break;
            
            /* wait until we have at least one diff-variance */
            if (sigma2Max == 0.0) continue;

            double epsAbs = relTol * bestAbs;                    // ε relative
            int    k      = inst.getItems();
            double logS   = k * Math.log(2.0);                   // |S| = 2^k
            double gamma  = (epsAbs-delta)*(epsAbs-delta)
                          / (2.0 * sigma2Max);
            int    Nld    = (int)Math.ceil( (logS - Math.log(alpha)) / gamma );

            System.out.printf("DEBUG LD:  N=%d σ̂²max=%.3g γ̂=%.3g N_LD=%d%n",
                              N, sigma2Max, gamma, Nld);

            if (N >= Nld) break;

            System.out.println("DEBUG LD: increasing N to " + (N<<1));
            N <<= 1;
            reset();
        }

        /* ---------- Phase 1 : CLT gap loop (both estimators) ---------- */
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
            double s2_M = varianceV(vBar);

            /*----- estimator 1 : fresh evaluation of incumbent ----------*/
            double[] eval = freshEvaluation(Nprime);
            double   gHat = eval[0];
            double   varN = eval[1] / Nprime;

            double point1 = vBar - gHat;
            double varTot = varN + s2_M/M;
            double z      = NormalDist.inverseF01(0.95);
            double gap1   = point1 + z*Math.sqrt(varTot);        // 95 % bound

            /*----- estimator 2 : uses small-scenario samples ------------*/
            double gBarN = 0.0;
            for (int m=0;m<M;m++)
                gBarN += computeSampleAverage(
                           reps.get(bestIdx).optimalKnapsack,
                           smallScen.get(m));
            gBarN /= M;

            double diffMean = vBar - gBarN;          // point estimator 2

            double s2bar = 0.0;
            for (int m=0;m<M;m++) {
                double gm = computeSampleAverage(
                               reps.get(bestIdx).optimalKnapsack,
                               smallScen.get(m));
                double vm = reps.get(m).milpSolutionValue;
                s2bar += Math.pow( (vm-gm) - diffMean , 2 );
            }
            s2bar = s2bar / (M*(M-1));               // \bar S²_M / M
            double gap2 = diffMean + z*Math.sqrt(s2bar);

            System.out.printf(
              "DEBUG GAP: rep=%2d  ĝ=%.6f  gap1=%.3e  gap2=%.3e  rel1=%.3e%n rel2=%.3e%n",
              M, gHat, gap1, gap2, gap1/gHat, gap2/gHat);

            if (gap1/gHat < relTol /*&& gap2/gHat < relTol*/) {
                finalGap1Rel = gap1/gHat;
                finalGap2Rel = gap2/gHat;
                finalEvalMean= gHat;
                break;
            }
            if (M >= Mmax) {
                finalGap1Rel = gap1/gHat;
                finalGap2Rel = gap2/gHat;
                finalEvalMean= gHat;
                System.out.println("DEBUG: reached Mmax");
                break;
            }
            try {
               makeReplication(N);
            } catch (IloException e) {
               // TODO Auto-generated catch block
               e.printStackTrace();
            }
            if (System.currentTimeMillis()-wall0 > wallMs)
                throw new RuntimeException("wall-clock limit");
        }

        long sec = (System.currentTimeMillis()-wall0)/1000;
        System.out.printf("TERMINATE: relGap1 %.3e  relGap2 %.3e  "
                          +"rep=%d  N=%d  time=%ds%n",
                          finalGap1Rel, finalGap2Rel, reps.size(), N, sec);

        /* -------- build solved-instance --------------------------------*/
        SKPGenericDistributionScenarioBasedSolvedInstance inc = reps.get(bestIdx);
        return new SKPGenericDistributionSAASolvedInstance(
                inst,
                inc.optimalKnapsack,
                finalEvalMean,
                sec*1000.0,
                finalGap1Rel,
                finalGap2Rel,
                N,
                Nprime,
                reps.size());
    }

    /*==================== helpers ====================================*/
    private void makeReplication(int N) throws IloException {

       SKPGenericDistributionScenarioBased saa =
              new SKPGenericDistributionScenarioBased(inst, N, rng);
       SKPGenericDistributionScenarioBasedSolvedInstance sol = saa.solve(Nprime);

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

       /* ------------- incumbent update test ---------------------------- */
       if (sol.simulatedSolutionValueMean1 > bestMean1) {
           bestMean1 = sol.simulatedSolutionValueMean1;
           bestAbs   = Math.abs(bestMean1);
           bestIdx   = reps.size() - 1;
           System.out.printf("DEBUG INC: new incumbent at rep %d "
                             +"Mean1=%.6f%n", reps.size(), bestMean1);

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
        double s=0;
        for (var r:reps)
            s += (r.milpSolutionValue-vBar)*(r.milpSolutionValue-vBar);
        return s/(reps.size()*(reps.size()-1));
    }

    private double computeSampleAverage(int[] knap, double[][] scen) {
        double val=0;
        for (int i=0;i<knap.length;i++)
            if (knap[i]==1) val += inst.getExpectedValues()[i];
        for (double[] row:scen) {
            double wt=0;
            for (int i=0;i<knap.length;i++)
                if (knap[i]==1) wt += row[i];
            val -= inst.getShortageCost()*Math.max(0, wt-inst.getCapacity());
        }
        return val/scen.length;
    }

    private void reset() {
        reps.clear();
        smallScen.clear();
        sigma2Max=0; sumV=0; bestIdx=-1;
        bestMean1=Double.NEGATIVE_INFINITY; bestAbs=0;
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
        return var / N;                     // unbiased N-denominator is fine here
    }


    /* ---------------- quick test ------------------------------------ */
    public static void main(String[] args) throws IloException {
        SKPGenericDistribution inst =
              SKPGenericDistribution.getTestInstanceLarge();

        SKPGenericDistributionSAA_LD solver =
              new SKPGenericDistributionSAA_LD(inst);

        SKPGenericDistributionSAASolvedInstance res = solver.solve();
        System.out.println("\nJSON RESULT:\n"+
                           GSONUtility.printInstanceAsJSON(res));
    }
}
