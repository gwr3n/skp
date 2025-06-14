package skp.saa;

import java.util.ArrayList;
import java.util.List;

import ilog.concert.IloException;
import skp.instance.SKPMultinormal;
import skp.milp.SKPMultinormalScenarioBased;
import skp.milp.instance.SKPMultinormalScenarioBasedSolvedInstance;
import skp.saa.instance.SKPMultinormalSAASolvedInstance;
import skp.sim.SimulateMultinormal;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

/**  SAA + large-deviation sample-size rule (Section 2) +
 *   two CLT gap estimators (Section 3) for MULTINORMAL weights.       */
public final class SKPMultinormalSAA_LD {

    /* ---------- user parameters ---------- */
    private static final double alpha  = 0.05;
    private static final double delta  = 0.0;
    private static final int    N0     = 64;
    private static final int    Mmax   = 100;
    private static final int    Nprime = 10_000;
    private static final int    warmUp = 32;
    private static final double relTol = 1e-2;
    private static final long   wallMs = 10 * 60_000L;

    /* ---------- members ---------- */
    private final SKPMultinormal inst;
    private final RandomStream   rng;

    private final List<SKPMultinormalScenarioBasedSolvedInstance> reps = new ArrayList<>();
    private final List<double[][]> smallScenarios = new ArrayList<>();   // store N scenarios

    private double sigma2Max = 0.0;
    private double sumV      = 0.0;
    private int    bestIdx   = -1;
    private double bestMean1 = Double.NEGATIVE_INFINITY;
    private double bestAbs   = 0.0;

    public SKPMultinormalSAA_LD(SKPMultinormal inst) {
        this.inst = inst;
        MRG32k3aL r = new MRG32k3aL();
        r.setSeed(new long[]{12345,24513,24531,42531,35124,32451});
        this.rng = r;
    }

    /* ================================================================= */
    public SKPMultinormalSAASolvedInstance solve() throws IloException {

        long wall0 = System.currentTimeMillis();
        int  N     = N0;

        /* --------- Phase 0: enlarge N until LD condition satisfied ------*/
        while (true) {
            makeReplication(N);

            if (reps.size() >= Mmax) break;
            
            /* wait until we have at least one diff-variance */
            if (sigma2Max == 0.0) continue;

            double epsAbs = relTol * bestAbs;              // ε = ρ·|ĝ|
            int    k      = inst.getItems();
            double logS   = k * Math.log(2.0);             // |S|=2^k
            double gamma  = (epsAbs - delta)*(epsAbs - delta) / (2.0*sigma2Max);
            int Nld       = (int)Math.ceil( (logS - Math.log(alpha)) / gamma );

            System.out.printf("DEBUG LD: N=%d σ̂²max=%.4g γ̂=%.4g N_LD=%d%n",
                              N, sigma2Max, gamma, Nld);

            if (N >= Nld) break;

            System.out.println("DEBUG LD: increasing N to " + (N<<1));
            N <<= 1;
            reset();
        }

        /* --------- Phase 1: CLT gap loop (gap1 & gap2) ------------------*/
        double gap1Rel = Double.POSITIVE_INFINITY, gap2Rel = Double.POSITIVE_INFINITY;
        double finalMean = Double.NaN;

        while (true) {

            if (reps.size() < warmUp) { makeReplication(N); continue; }

            int    M  = reps.size();
            double vBar = sumV / M;
            double s2_M = varianceV(vBar);

            /* ---- gap1 (fresh evaluation) ---- */
            double[] eval = freshEvaluation(Nprime);
            double   gHat = eval[0];
            double   s2Np = eval[1] / Nprime;

            double z   = NormalDist.inverseF01(0.95);
            double gap1= (vBar - gHat) + z*Math.sqrt(s2Np + s2_M/M);

            /* ---- gap2 (Section 3, alternative) ---- */
            double gBarN = 0.0;
            for (int m=0; m<M; m++)
                gBarN += sampleAverage(reps.get(bestIdx).optimalKnapsack,
                                       smallScenarios.get(m));
            gBarN /= M;

            double diff   = vBar - gBarN;          // point estimator
            double s2bar  = 0.0;
            for (int m=0; m<M; m++) {
                double gm = sampleAverage(reps.get(bestIdx).optimalKnapsack,
                                          smallScenarios.get(m));
                double vm = reps.get(m).milpSolutionValue;
                s2bar += Math.pow( (vm - gm) - diff , 2 );
            }
            s2bar /= (M * (M-1));                  //  \bar S²_M / M(M−1)
            double gap2 = diff + z*Math.sqrt(s2bar);

            System.out.printf("DEBUG GAP: rep=%2d ĝ=%.6f "
                              +"gap1=%.3e gap2=%.3e rel1=%.3e%n rel2=%.3e%n",
                              M, gHat, gap1, gap2, gap1/gHat, gap2 / gHat);

            gap1Rel = gap1 / gHat;
            gap2Rel = gap2 / gHat;
            finalMean = gHat;

            if (gap1Rel < relTol /*&& gap2Rel < relTol*/) break;
            if (M >= Mmax) {
                System.out.println("DEBUG: reached Mmax");
                break;
            }
            makeReplication(N);
            if (System.currentTimeMillis()-wall0 > wallMs)
                throw new RuntimeException("wall-clock limit");
        }

        long sec = (System.currentTimeMillis()-wall0)/1000;
        System.out.printf("TERMINATE: relGap1 %.3e  relGap2 %.3e  "
                          +"rep=%d N=%d time=%ds%n",
                          gap1Rel, gap2Rel, reps.size(), N, sec);

        /* ------------ build solved-instance --------------------------- */
        SKPMultinormalScenarioBasedSolvedInstance inc = reps.get(bestIdx);
        return new SKPMultinormalSAASolvedInstance(
                inst,
                inc.optimalKnapsack,
                finalMean,
                sec*1000.0,
                gap1Rel,
                gap2Rel,
                N,
                Nprime,
                reps.size());
    }

    /* ==================== helpers ====================================*/
    private void makeReplication(int N) throws IloException {

        SKPMultinormalScenarioBased saa =
                new SKPMultinormalScenarioBased(inst, N, rng);
        SKPMultinormalScenarioBasedSolvedInstance sol =
                saa.solve(Nprime);

        reps.add(sol);
        smallScenarios.add(saa.scenarios);          // store N-scenario matrix

        sumV      += sol.milpSolutionValue;
        sigma2Max  = Math.max(sigma2Max, sol.simulatedSolutionValueVariance1);

        if (sol.simulatedSolutionValueMean1 > bestMean1) {
            bestMean1 = sol.simulatedSolutionValueMean1;
            bestAbs   = Math.abs(bestMean1);
            bestIdx   = reps.size()-1;
            System.out.printf("DEBUG INC: new incumbent at rep %d "
                              +"Mean1=%.6f%n", reps.size(), bestMean1);
        }
    }

    private double[] freshEvaluation(int nSim) {
        SimulateMultinormal sim = new SimulateMultinormal(inst);
        return sim.simulateMeanVariance(reps.get(bestIdx).optimalKnapsack,
                                        nSim, rng);
    }

    private double varianceV(double vBar) {
        double s = 0.0;
        for (var r: reps)
            s += (r.milpSolutionValue - vBar)
                 * (r.milpSolutionValue - vBar);
        return s / (reps.size() * (reps.size()-1));
    }

    private double sampleAverage(int[] knap, double[][] scen) {
        double value = 0.0;
        for (int i=0;i<knap.length;i++)
            if (knap[i]==1) value += inst.getExpectedValues()[i];

        for (double[] row: scen) {
            double wt=0;
            for (int i=0;i<knap.length;i++)
                if (knap[i]==1) wt += row[i];
            value -= inst.getShortageCost()*Math.max(0, wt-inst.getCapacity());
        }
        return value / scen.length;
    }

    private void reset() {
        reps.clear();
        smallScenarios.clear();
        sigma2Max=0; sumV=0; bestIdx=-1;
        bestMean1=Double.NEGATIVE_INFINITY; bestAbs=0;
        rng.resetStartStream();
    }

    /* ---------------- simple test driver ----------------------------- */
    public static void main(String[] args) throws IloException {

        SKPMultinormal inst = SKPMultinormal.getTestInstance();
        SKPMultinormalSAA_LD solver = new SKPMultinormalSAA_LD(inst);

        SKPMultinormalSAASolvedInstance res = solver.solve();
        System.out.println("\nJSON RESULT:\n"
                           + GSONUtility.printInstanceAsJSON(res));
    }
}
