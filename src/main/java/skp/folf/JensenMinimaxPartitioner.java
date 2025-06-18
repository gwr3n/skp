/******************************************************************
 * JensenMinimaxPartitioner
 * --------------------------------------------------------------
 * Builds an optimal (minimax) Jensen lower bound with a given
 * number of segments for the complementary first–order loss
 * function of a standard Normal random variable Z~N(0,1).
 *
 * References:
 *   R. Rossi, S. A. Tarim, B. Hnich, S. Prestwich,
 *   "Piecewise linear lower and upper bounds for the standard normal first order loss function",
 *   Applied Mathematics & Computation 231 (2014) 489-502 — Sec. 4.2.1
 *
 * Requires:
 *   – SSJ  (umontreal.ssj)
 *   – Apache Commons-Math 3.x
 ******************************************************************/

package skp.folf;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import umontreal.ssj.probdist.NormalDist;

/**
 * Computes the minimax Jensen lower bound for the complementary first-order loss function
 * of a standard normal random variable, using a specified number of segments.
 * The algorithm follows the approach in Rossi et al. (2014), Sec. 4.2.1.
 */
public final class JensenMinimaxPartitioner {

    /**
     * Result container for the minimax partitioning.
     * Holds the segment probabilities, conditional means, and the maximum error.
     */
    public static final class Result {
        /** Segment probabilities \( p_i \) */
        public final double[] p;
        /** Conditional means \( E[Z|\Omega_i] \) */
        public final double[] expect;
        /** Maximum error of the Jensen lower bound */
        public final double error;

        private Result(double[] p, double[] e, double err) {
            this.p = p;
            this.expect = e;
            this.error = err;
        }
    }

    /**
     * Converts a vector of log-probabilities to breakpoints (cumulative log-sums).
     * This is used to ensure positivity and monotonicity of the breakpoints.
     * See Rossi et al. (2014), Eq. (19).
     */
    private static double[] toBreakpoints(double[] y) {
        double[] bp = new double[y.length];
        double cum = Math.log(0); // Start from log(0) == negative infinity
        for (int i = 0; i < y.length; ++i) {
            cum = logSumExp(cum, y[i]); // Numerically stable log-sum
            bp[i] = cum;
        }
        return bp;
    }

    /**
     * Numerically stable computation of log(exp(loga) + exp(logb)).
     * Prevents underflow/overflow in log-space.
     */
    private static double logSumExp(double loga, double logb) {
        if (Double.isInfinite(loga)) return logb;
        if (Double.isInfinite(logb)) return loga;
        double max = Math.max(loga, logb);
        return max + Math.log(Math.exp(loga - max) + Math.exp(logb - max));
    }

    /**
     * Computes the minimax Jensen lower bound partition for a given number of segments.
     * @param segments Number of segments (must be >= 2)
     * @return Result object containing probabilities, means, and error
     */
    public static Result compute(int segments) {
        if (segments < 2)
            throw new IllegalArgumentException("segments must be ≥ 2");

        final int N = segments - 1; // Number of intervals

        // Special case: N=1 (2 segments), trivial partition
        if (N == 1) {
            double[] p = { 1.0 };
            double[] Ei = { 0.0 };
            double err = 1.0 / Math.sqrt(2.0 * Math.PI); // Standard normal density at 0
            return new Result(p, Ei, err);
        }

        // Special case: N=2 (3 segments), explicit solution
        if (N == 2) {
            double[] p = { 0.5, 0.5 };
            double phi0 = NormalDist.density01(0.0);
            double Ei = -phi0 / 0.5;
            double[] E = { Ei, -Ei };

            double PhiEi = NormalDist.cdf01(-Ei);
            double loss = NormalDist.density01(Ei) + (-Ei) * PhiEi;

            double lbVal = 0.5 * (-Ei) + phi0;

            double err = loss - lbVal;
            return new Result(p, E, err);
        }

        // For N > 2, use optimization as in Rossi et al. (2014), Sec. 4.2.1
        final int mPos = (N - 1) / 2; // Number of positive breakpoints
        final boolean addZero = (N % 2 == 0); // Even N: add zero breakpoint
        final int dim = mPos; // Optimization dimension

        // Objective function: minimize the squared difference of errors (see Eq. (21))
        Objective fObj = new Objective(N, mPos, addZero);

        // Initial guess: log(0.5 + j) for each positive breakpoint
        double[] y0 = new double[dim];
        for (int j = 0; j < dim; ++j)
            y0[j] = Math.log(0.5 + j);

        // Set up the Nelder-Mead simplex optimizer (derivative-free)
        SimplexOptimizer opt = new SimplexOptimizer(
                new SimplePointChecker<>(1e-12, 1e-12));

        // Run the optimization to find the best breakpoints
        PointValuePair sol = opt.optimize(
                new MaxIter(100000), new MaxEval(100000),
                new ObjectiveFunction(fObj),
                GoalType.MINIMIZE,
                new InitialGuess(y0),
                new NelderMeadSimplex(dim)
        );

        // Convert optimized log-probabilities to breakpoints
        double[] bpPos = toBreakpoints(sol.getPoint());

        // Build the full set of breakpoints (symmetric about zero)
        int B = 2 * mPos + (addZero ? 1 : 0);
        double[] bInt = new double[B];
        for (int j = 0; j < mPos; ++j)
           bInt[j] = -bpPos[mPos - 1 - j];
        if (addZero)
           bInt[mPos] = 0.0;
        for (int j = 0; j < mPos; ++j)
            bInt[mPos + (addZero ? 1 : 0) + j] = bpPos[j];

        // Compute interval bounds a[i], b[i] for each segment
        double[] a = new double[N], b = new double[N];
        a[0] = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < N - 1; ++i) {
            b[i] = bInt[i];
            a[i + 1] = bInt[i];
        }
        b[N - 1] = Double.POSITIVE_INFINITY;

        // Compute probabilities, conditional means, and error for each segment
        double[] p = new double[N], Ei = new double[N];
        double err = 0.0;
        for (int i = 0; i < N; ++i) {
            double Phi_ai = (i == 0) ? 0.0 : NormalDist.cdf01(a[i]);
            double Phi_bi = (i == N - 1) ? 1.0 : NormalDist.cdf01(b[i]);
            double phi_ai = (i == 0) ? 0.0 : NormalDist.density01(a[i]);
            double phi_bi = (i == N - 1) ? 0.0 : NormalDist.density01(b[i]);

            p[i] = Phi_bi - Phi_ai; // Probability of segment i
            Ei[i] = (phi_ai - phi_bi) / p[i]; // Conditional mean in segment i

            // Compute the error for this segment (see Eq. (20))
            double loss = NormalDist.density01(Ei[i]) + Ei[i] * NormalDist.cdf01(Ei[i]);
            double lb = Phi_bi * Ei[i] + phi_bi;
            double ei = loss - lb;
            if (ei > err) err = ei; // Track the maximum error
        }
        return new Result(p, Ei, err);
    }

    /**
     * Objective function for the minimax optimization.
     * The goal is to minimize the squared difference of errors across segments,
     * enforcing the minimax property (see Rossi et al. (2014), Eq. (21)).
     */
    private static final class Objective implements MultivariateFunction {
        private final int N, mPos;
        private final boolean addZero;

        Objective(int N, int mPos, boolean addZero) {
            this.N = N;
            this.mPos = mPos;
            this.addZero = addZero;
        }

        @Override
        public double value(double[] y) {
            // Convert log-probabilities to breakpoints
            double[] bpPos = toBreakpoints(y);
            int B = 2 * mPos + (addZero ? 1 : 0);
            double[] bInt = new double[B];
            for (int j = 0; j < mPos; ++j)
               bInt[j] = -bpPos[mPos - 1 - j];
            if (addZero)
               bInt[mPos] = 0.0;
            for (int j = 0; j < mPos; ++j)
                bInt[mPos + (addZero ? 1 : 0) + j] = bpPos[j];

            // Compute interval bounds
            double[] a = new double[N], b = new double[N], e = new double[N];
            a[0] = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < N - 1; ++i) {
                b[i] = bInt[i];
                a[i + 1] = bInt[i];
            }
            b[N - 1] = Double.POSITIVE_INFINITY;

            // Compute error for each segment
            for (int i = 0; i < N; ++i) {
                double Phi_ai = (i == 0) ? 0 : NormalDist.cdf01(a[i]);
                double Phi_bi = (i == N - 1) ? 1 : NormalDist.cdf01(b[i]);
                double phi_ai = (i == 0) ? 0 : NormalDist.density01(a[i]);
                double phi_bi = (i == N - 1) ? 0 : NormalDist.density01(b[i]);
                double p = Phi_bi - Phi_ai;
                double E = (phi_ai - phi_bi) / p;
                double loss = NormalDist.density01(E) + E * NormalDist.cdf01(E);
                double lb = Phi_bi * E + phi_bi;
                e[i] = loss - lb;
            }
            // Objective: sum of squared differences from the first segment's error
            double obj = 0;
            for (int i = 1; i < N; ++i) {
                double d = e[i] - e[0];
                obj += d * d;
            }
            return obj;
        }
    }

    /**
     * Self-test routine: checks computed results against reference values from Rossi et al. (2014).
     * Prints results and verifies accuracy.
     */
    public static void selfTest() {
       final double[] REF_ERRORS = getErrors();
       final double EPS = 1e-9;
       System.out.println("----- Jensen minimax self-test -----");

       for (int W = 1; W <= 10; ++W) {
           Result r = compute(W + 1);
           double[] expectedMeans = getMeans(W);
           double[] expectedProbs = getProbabilities(W);

           double errorDiff = Math.abs(r.error - REF_ERRORS[W - 1]);
           boolean errorCheck = (errorDiff <= EPS);

           boolean meansCheck = true;
           for (int i = 0; i < expectedMeans.length; i++) {
               if (Math.abs(r.expect[i] - expectedMeans[i]) > EPS) {
                   meansCheck = false;
                   break;
               }
           }

           boolean probsCheck = true;
           for (int i = 0; i < expectedProbs.length; i++) {
               if (Math.abs(r.p[i] - expectedProbs[i]) > EPS) {
                   probsCheck = false;
                   break;
               }
           }

           System.out.printf(
               "W=%2d  error=%.15f  ref=%.15f  diff=%8.2e %s  meansCheck=%s  probsCheck=%s%n",
               W, r.error, REF_ERRORS[W - 1], errorDiff,
               errorCheck ? "OK" : "FAIL",
               meansCheck ? "OK" : "FAIL",
               probsCheck ? "OK" : "FAIL"
           );

           if (!meansCheck || !probsCheck) {
               System.err.printf("Test case failed for W = %d:\n", W);
               if (!meansCheck) {
                   System.err.println("  Expected means: " + arrayToString(expectedMeans));
                   System.err.println("  Computed means: " + arrayToString(r.expect));
               }
               if (!probsCheck) {
                   System.err.println("  Expected probabilities: " + arrayToString(expectedProbs));
                   System.err.println("  Computed probabilities: " + arrayToString(r.p));
               }
               System.exit(1);
           }
       }

       System.out.println("All tests passed.");
   }

   /**
    * Utility: formats a double array as a string for output.
    */
   private static String arrayToString(double[] array) {
       StringBuilder sb = new StringBuilder();
       sb.append("[");
       for (int i = 0; i < array.length; i++) {
           sb.append(String.format("%.8f", array[i]));
           if (i < array.length - 1) {
               sb.append(", ");
           }
       }
       sb.append("]");
       return sb.toString();
   }

   //----------------------------------------------------------------------------
   // Reference Getters
   //----------------------------------------------------------------------------
   
   /**
    * Returns reference errors for the given number of partitions,
    * as reported in Rossi et al. (2014), Table 2.
    */
   public static double[] getErrors() {
      return new double[] {
            0.3989422804014327,
            0.1206560496714961,
            0.05784405029198253,
            0.033905164962384104,
            0.022270929512393414,
            0.01574607463566398,
            0.011721769576577057,
            0.00906528789647753,
            0.007219916411227892,
            0.005885974956458359
        };
   }
   
   /**
    * Returns reference conditional means for the given number of partitions,
    * as reported in Rossi et al. (2014), Table 2.
    */
   public static double[] getMeans(int partitions) {
       switch (partitions) {
           case 1:
               return new double[] {0};
           case 2:
               return new double[] {-0.7978845608028654, 0.7978845608028654};
           case 3:
               return new double[] {-1.1850544278068644, 0, 1.1850544278068644};
           case 4:
               return new double[] {-1.4353532729205845, -0.41522324304905966, 0.41522324304905966, 1.4353532729205845};
           case 5:
               return new double[] {-1.6180463502161044, -0.6914240068499904, 0, 0.6914240068499903, 1.6180463502161053};
           case 6:
               return new double[] {-1.7608020666235031, -0.8960107374480083, -0.28188851144130117, 0.28188851144130117, 0.8960107374480083, 1.7608020666235031};
           case 7:
               return new double[] {-1.8773528492652836, -1.0572304450884658, -0.4934048390251067, 0, 0.4934048390251065, 1.0572304450884658, 1.8773528492652836};
           case 8:
               return new double[] {-1.9754694729585056, -1.1895340795157716, -0.6615516528578579, -0.213586638906901, 0.213586638906901, 0.6615516528578579, 1.1895340795157716, 1.9754694729585056};
           case 9:
               return new double[] {-2.059957433491476, -1.30127090280595, -0.8004000560466271, -0.3845969617811554, 0, 0.3845969617811554, 0.8004000560466271, 1.30127090280595, 2.059957433491476};
           case 10:
               return new double[] {-2.133986195498256, -1.3976822972668839, -0.918199946431143, -0.5265753462727588, -0.17199013069262026, 0.17199013069262026, 0.5265753462727588, 0.918199946431143, 1.3976822972668839, 2.133986195498256};
           default:
               throw new IllegalArgumentException("Segments not supported in test.");
       }
   }

   /**
    * Returns reference probabilities for the given number of partitions,
    * as reported in Rossi et al. (2014), Table 2.
    */
   public static double[] getProbabilities(int partitions) {
       switch (partitions) {
           case 1:
               return new double[] {1};
           case 2:
               return new double[] {0.5, 0.5};
           case 3:
               return new double[] {0.28783338731597996, 0.4243332253680401, 0.28783338731597996};
           case 4:
               return new double[] {0.18755516774758485, 0.31244483225241515, 0.31244483225241515, 0.18755516774758485};
           case 5:
               return new double[] {0.1324110437406592, 0.23491250409192982, 0.26535290433482195, 0.23491250409192987, 0.13241104374065915};
           case 6:
               return new double[] {0.09877694599482933, 0.18223645973091096, 0.2189865942742597, 0.2189865942742597, 0.18223645973091096, 0.09877694599482933};
           case 7:
               return new double[] {0.07669891602586965, 0.14538182479573014, 0.18144834397745296, 0.19294183040189444, 0.18144834397745302, 0.14538182479573014, 0.07669891602586965};
           case 8:
               return new double[] {0.061394553470121016, 0.11872108750901467, 0.15205073490726895, 0.16783362411359537, 0.16783362411359537, 0.15205073490726895, 0.11872108750901467, 0.061394553470121016};
           case 9:
               return new double[] {0.05033061540430428, 0.09884442068481658, 0.12900389832341652, 0.1460368860056377, 0.15156835916364986, 0.1460368860056377, 0.12900389832341652, 0.09884442068481658, 0.05033061540430428};
           case 10:
               return new double[] {0.04206108420763477, 0.0836356495308449, 0.11074334596058821, 0.1276821455299152, 0.13587777477101692, 0.13587777477101692, 0.1276821455299152, 0.11074334596058821, 0.0836356495308449, 0.04206108420763477};
           default:
               throw new IllegalArgumentException("Segments not supported in test.");
       }
   }

    /**
     * Main entry point. If no arguments, runs self-test.
     * If an integer argument is given, computes and prints the optimal Jensen lower bound
     * for that number of segments.
     */
    public static void main(String[] args) {
        if (args.length == 0) {
            selfTest();
            return;
        }
        int segments = Integer.parseInt(args[0]);
        Result r = compute(segments);
        System.out.printf("%nOptimal Jensen lower bound with %d segments%n", segments);
        System.out.println(" i       p_i        E[Z|Ω_i]");
        for(int i = 0; i < r.p.length; ++i)
            System.out.printf("%2d  %11.8f  %11.8f%n", i+1, r.p[i], r.expect[i]);
        System.out.printf("Common max error e_W = %.15f%n", r.error);
    }
}
