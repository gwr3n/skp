package skp.folf;

public class LinearisationFactory {

   public static int[] chooseLinearisationParameters(final double epsilon, double Vmax, double c)
   {
      JensenPartitioner partitioner;
      switch(PiecewiseStandardNormalFirstOrderLossFunction.getPartitioner()) {
      case PARTITIONER.UNIFORM:
         partitioner = new JensenUniformPartitioner();
         break;
      case PARTITIONER.MINIMAX:
      default:
         partitioner = new JensenMinimaxPartitioner();
         break;
      }

      if (epsilon <= 0.0)
         throw new IllegalArgumentException("epsilon must be positive");

      /* ---------- constants that depend only on the instance ---------- */
      double Smax = Math.sqrt(Vmax);

      /* ================================================================
       * 1.  Choose W  (first-order-loss approximation)
       * ================================================================ */
      final double rhsLoss = epsilon / (2.0 * c * Smax);   // ε /(2 c S_max)

      int    W = 1;                        // current number of partitions
      Result part = partitioner.compute(W + 1);

      /* ---- exponential search to find an upper bound that satisfies condition */
      while (part.error > rhsLoss) {
         W <<= 1;                                        // double W
         part = partitioner.compute(W + 1);
         if (W > (1 << 22))                              // crude overflow guard
            throw new IllegalStateException("W exploded (> 4,000,000)");
      }

      /* ---- binary search for the *minimal* feasible W ---------------- */
      int low = W >> 1, high = W;
         while (low + 1 < high) {
            int mid = (low + high) >>> 1;                   // midpoint
         Result midPart = partitioner.compute(mid + 1);
         if (midPart.error <= rhsLoss) {
            high  = mid;
            part  = midPart;
         } else {
            low = mid;
         }
         }
         W = high;
         final double dotErr = part.error;                   // ẑe_W just chosen
         final int    wSegments = W + 1;

         /* ---------- compute  A_max  for this partition ------------------ */
         final double[] p  = part.p;     // length W
         final double[] mu = part.expect;     // length W
         double Amax = 0.0;
         for (int i = 0; i < W; ++i) {
            double sum = 0.0;
            for (int k = i; k < W; ++k)
               sum += p[k] * mu[k];
            if (sum > Amax) Amax = sum;
         }

         /* ================================================================
          * 2.  Choose Q  (square-root approximation)
          * ================================================================ */
         /* RHS of  δ_Q ≤ ε /(2 c (A_max + dotErr)) */
         double rhsDelta = epsilon / (2.0 * c * (Amax + dotErr));

         /* δ_Q = √(Vmax / Q)/4  ⇒  Q ≥ Vmax /(16 · rhsDelta²) */
         double qMin = Vmax / (16.0 * rhsDelta * rhsDelta);
         int Q = (int)Math.ceil(qMin);
         if (Q < 1) Q = 1;                                   // sanity floor

         /* Optional run-time assertion (disable in production) */
         assert Math.sqrt(Vmax / Q) / 4.0 <= rhsDelta + 1e-12;

         //System.out.println("Chosen linearisation parameters: W = " + W + "; s = " + (Vmax/Q));

         return new int[] { wSegments, Q };
   }
}
