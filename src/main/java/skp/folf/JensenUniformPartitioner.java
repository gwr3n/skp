// language: java
package skp.folf;

import umontreal.ssj.probdist.NormalDist;

public final class JensenUniformPartitioner extends JensenPartitioner {


   // --- Numerical helpers ---

   private static final double LOG_SQRT_2PI = Math.log(Math.sqrt(2.0 * Math.PI));

   private static double phi(double z) {
      return NormalDist.density01(z);
   }

   private static double Phi(double z) {
      return NormalDist.cdf01(z);
   }

   // Complementary tail via symmetry: barΦ(z) = Φ(-z)
   private static double barPhi(double z) {
      return NormalDist.cdf01(-z);
   }

   // log φ(z) = -0.5 z^2 - log(√(2π))
   private static double logPhi(double z) {
      return -0.5 * z * z - LOG_SQRT_2PI;
   }

   // Stable difference Φ(b) - Φ(a); uses symmetry when both are in right tail.
   // Falls back to Simpson for extremely small intervals.
   private static double cdfDiffStable(double a, double b) {
      if (a == Double.NEGATIVE_INFINITY) return (b == Double.POSITIVE_INFINITY) ? 1.0 : Phi(b);
      if (b == Double.POSITIVE_INFINITY) return (a == Double.NEGATIVE_INFINITY) ? 1.0 : barPhi(a);
      if (a >= 0.0 && b >= 0.0) {
         // right tail: use Φ(-a) - Φ(-b)
         double d = Phi(-a) - Phi(-b);
         if (d <= 0.0) d = 0.0; // guard tiny negative from rounding
         // Fallback if cancellation still bites:
         if (d < 1e-16) {
            double m = 0.5 * (a + b);
            double dSimpson = (b - a) / 6.0 * (phi(a) + 4.0 * phi(m) + phi(b));
            return dSimpson;
         }
         return d;
      } else {
         // general case: direct difference is safe (not both near 1)
         double d = Phi(b) - Phi(a);
         if (d <= 0.0) {
            // Fallback Simpson if rounding inverted ordering
            double m = 0.5 * (a + b);
            return Math.max(0.0, (b - a) / 6.0 * (phi(a) + 4.0 * phi(m) + phi(b)));
         }
         return d;
      }
   }

   // Stable difference φ(a) - φ(b) using logs + expm1
   private static double densityDiffStable(double a, double b) {
      if (a == Double.NEGATIVE_INFINITY) return -phi(b); // φ(-∞)=0
      if (b == Double.POSITIVE_INFINITY) return phi(a);  // φ(+∞)=0
      double la = logPhi(a);
      double lb = logPhi(b);
      if (la >= lb) {
         double delta = lb - la; // <= 0
         // exp(la) * (1 - exp(delta)) = -exp(la) * expm1(delta)
         return -Math.exp(la) * Math.expm1(delta);
      } else {
         double delta = la - lb; // < 0
         // -( exp(lb) * (1 - exp(delta)) ) = Math.exp(lb) * (-expm1(delta)) with sign
         return Math.exp(lb) * Math.expm1(delta);
      }
   }

   // Mills ratio for right tail: r(a) = φ(a) / Φ(-a)
   private static double millsRight(double a) {
      double numLog = logPhi(a);
      double den = barPhi(a);        // = Φ(-a), small but stable
      return Math.exp(numLog) / den; // safe even for very small num; ratio is moderate
   }

   // Mills ratio for left tail: λ(b) = φ(b) / Φ(b)
   private static double millsLeft(double b) {
      double numLog = logPhi(b);
      double den = Phi(b);           // small but stable for b << 0
      return Math.exp(numLog) / den;
   }

   // Truncated mean E[Z | a < Z < b] with stable tails and mid-interval case
   private static double truncatedMean(double a, double b, double p) {
      if (a == Double.NEGATIVE_INFINITY) {
         if (b == Double.POSITIVE_INFINITY) return 0.0; // whole normal
         // left tail mean: -φ(b)/Φ(b)
         return -millsLeft(b);
      }
      if (b == Double.POSITIVE_INFINITY) {
         // right tail mean: φ(a)/Φ(-a)
         return millsRight(a);
      }
      // mid-interval: (φ(a) - φ(b)) / p using stable numerator/denominator
      double num = densityDiffStable(a, b);
      return num / p;
   }


   /** Compute uniform partition */

   public Result compute(int segments) {
      if (segments < 2)
         throw new IllegalArgumentException("segments must be ≥ 2");

      final int N = segments - 1;

      // N == 1: single segment covering whole normal
      if (N == 1) {
         double[] p = { 1.0 };
         double[] Ei = { 0.0 };
         double err = 1.0 / Math.sqrt(2.0 * Math.PI); // unchanged
         return new Result(p, Ei, err);
      }

      // Uniform-probability breakpoints (quantiles)
      double[] breakpoints = new double[N - 1];
      for (int i = 1; i < N; ++i) {
         breakpoints[i - 1] = NormalDist.inverseF01((double) i / N);
      }

      // Segment bounds
      double[] a = new double[N];
      double[] b = new double[N];
      a[0] = Double.NEGATIVE_INFINITY;
      for (int i = 0; i < N - 1; ++i) {
         b[i] = breakpoints[i];
         a[i + 1] = breakpoints[i];
      }
      b[N - 1] = Double.POSITIVE_INFINITY;

      // Outputs
      double[] p = new double[N];
      double[] Ei = new double[N];

      // Compute probabilities and conditional means with stability
      double sumP = 0.0;
      for (int i = 0; i < N; ++i) {
         double Pi;  // Φ(b_i)
         if (i == 0) {
            // first segment: p = Φ(b)
            Pi = (N == 2) ? Phi(b[i]) : Phi(b[i]);
            p[i] = Pi;
         } else if (i == N - 1) {
            // last segment: p = 1 - Φ(a) = Φ(-a)
            Pi = 1.0;
            p[i] = barPhi(a[i]); // stable via symmetry
         } else {
            // mid segments: stable difference
            p[i] = cdfDiffStable(a[i], b[i]);
         }

         // Guard against extremely tiny or zero p due to rounding
         if (p[i] < 0.0) p[i] = 0.0;
         sumP += p[i];

         // Conditional mean
         Ei[i] = truncatedMean(a[i], b[i], p[i]);
      }

      // Renormalize probabilities to sum exactly 1 (important for large N)
      if (sumP > 0.0 && Math.abs(sumP - 1.0) > 1e-15) {
         double scale = 1.0 / sumP;
         for (int i = 0; i < N; ++i) p[i] *= scale;
      }



      // --- Error computation (unchanged formula, now with stable inputs) ---
      double err = 0.0;
      for (int i = 0; i < N; ++i) {
        double Phi_bi = (i == N - 1) ? 1.0 : Phi(b[i]);        
        double phi_bi = (i == N - 1) ? 0.0 : phi(b[i]);

        double loss = phi(Ei[i]) + Ei[i] * Phi(Ei[i]);
        double lb = Phi_bi * Ei[i] + phi_bi;
        double ei = loss - lb;
        if (ei > err) err = ei;
      }



      return new Result(p, Ei, err);
   }
   

   public static void selfTest() {
       System.out.println("----- Jensen uniform partition diagnostic -----");
       final int[] TEST_W = {2, 4, 8, 16, 32, 64};
       final double EPS = 1e-12;
       
       JensenUniformPartitioner partitioner = new JensenUniformPartitioner();
   
       for (int W : TEST_W) {
           Result r = partitioner.compute(W + 1);
           double sumP = 0.0;
           for (double pi : r.p) sumP += pi;
   
           boolean normOK = Math.abs(sumP - 1.0) < EPS;
           boolean monoOK = true;
           for (int i = 1; i < r.expect.length; ++i) {
               if (r.expect[i] <= r.expect[i - 1]) {
                   monoOK = false;
                   break;
               }
           }
   
           boolean symOK = true;
           for (int i = 0; i < r.p.length / 2; ++i) {
               if (Math.abs(r.p[i] - r.p[r.p.length - 1 - i]) > 1e-6 ||
                   Math.abs(r.expect[i] + r.expect[r.expect.length - 1 - i]) > 1e-6) {
                   symOK = false;
                   break;
               }
           }
           
           JensenMinimaxPartitioner minimaxPartitioner = new JensenMinimaxPartitioner();
           Result rMinimax = minimaxPartitioner.compute(W + 1);
           System.out.printf("W=%3d | sumP=%.2f norm=%s mono=%s sym=%s err=%.15f err_minimax=%.15f%n",
               W, sumP, normOK ? "OK" : "FAIL", monoOK ? "OK" : "FAIL", symOK ? "OK" : "FAIL", r.error, rMinimax.error);
   
           if (!normOK || !monoOK || !symOK) {
               System.err.println("Diagnostic failed for W=" + W);
               System.exit(1);
           }
       }
       System.out.println("All diagnostics passed.");
   }

   public static void main(String[] args) {
      if (args.length == 0) {
         selfTest();
         return;
      }
      int segments = Integer.parseInt(args[0]);
      JensenUniformPartitioner partitioner = new JensenUniformPartitioner();
      Result r = partitioner.compute(segments);
      System.out.printf("%nUniform Jensen lower bound with %d segments%n", segments);
      System.out.println(" i       p_i        E[Z|Ω_i]");
      for (int i = 0; i < r.p.length; ++i)
         System.out.printf("%2d  %11.8f  %11.8f%n", i + 1, r.p[i], r.expect[i]);
      System.out.printf("Common max error e_W = %.15f%n", r.error);
   }
}
