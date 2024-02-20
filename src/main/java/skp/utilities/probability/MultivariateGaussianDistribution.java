package skp.utilities.probability;

import java.util.Random;

import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import umontreal.ssj.probdist.NormalDist;

public class MultivariateGaussianDistribution {
   
   
   /** The mean vector. */
   public final double[] mu;
   /** The covariance matrix. */
   public final RealMatrix sigma;
   /** True if the covariance matrix is diagonal. */
   public final boolean diagonal;
   
   /** The dimension. */
   private int dim;
   /** The Cholesky decomposition of covariance matrix. */
   private RealMatrix sigmaL;
   
   /** The number of parameters. */
   private final int length;
   
   /**
    * Constructor.
    *
    * @param mean mean vector.
    * @param cov covariance matrix.
    */
   public MultivariateGaussianDistribution(double[] mean, RealMatrix cov) {
       if (mean.length != cov.getRowDimension()) {
           throw new IllegalArgumentException("Mean vector and covariance matrix have different dimension");
       }

       mu = mean;
       sigma = cov;

       diagonal = false;
       length = mu.length + mu.length * (mu.length + 1) / 2;

       init();
   }
   
   /**
    * Initialize the object.
    */
   private void init() {
       dim = mu.length;
       
       CholeskyDecomposition cholesky = new CholeskyDecomposition(sigma);
       sigmaL = cholesky.getL();
   }
   
   public int length() {
       return length;
   }

   public double[] mean() {
       return mu;
   }

   public RealMatrix cov() {
       return sigma;
   }
   
   /**
    * Algorithm from Alan Genz (1992) Numerical Computation of 
    * Multivariate Normal Probabilities, Journal of Computational and 
    * Graphical Statistics, pp. 141-149.
    *
    * The difference between returned value and the true value of the
    * CDF is less than 0.001 in 99.9% time. The maximum number of iterations
    * is set to 10000.
    */
   public double cdf(double[] x) {
       Random rnd = new Random(123456);
      
       if (x.length != dim) {
           throw new IllegalArgumentException("Sample has different dimension.");
       }

       int Nmax = 10000;
       double alph = NormalDist.inverseF01(0.999);
       double errMax = 0.001;

       double[] v = x.clone();
       sub(v, mu);

       double p = 0.0;
       double varSum = 0.0;

       // d is always zero
       double[] e = new double[dim];
       double[] f = new double[dim];
       e[0] = NormalDist.cdf01(v[0] / sigmaL.getEntry(0, 0));
       f[0] = e[0];

       double[] y = new double[dim];

       double err = 2 * errMax;
       int N;
       for (N = 1; err > errMax && N <= Nmax; N++) {
           double[] w = rnd.doubles(dim - 1).toArray();
           for (int i = 1; i < dim; i++) {
               y[i - 1] = NormalDist.inverseF01(w[i - 1] * e[i - 1]);
               double q = 0.0;
               for (int j = 0; j < i; j++) {
                   q += sigmaL.getEntry(i, j) * y[j];
               }

               e[i] = NormalDist.cdf01((v[i] - q) / sigmaL.getEntry(i, i));
               f[i] = e[i] * f[i - 1];
           }

           double del = (f[dim - 1] - p) / N;
           p += del;
           varSum = (N - 2) * varSum / N + del * del;
           err = alph * Math.sqrt(varSum);
       }

       return p;
   }
   
   /**
    * Element-wise subtraction of two arrays y = y - x.
    * @param y the minuend array.
    * @param x the subtrahend array.
    */
   public static void sub(double[] y, double[] x) {
       if (x.length != y.length) {
           throw new IllegalArgumentException(String.format("Arrays have different length: x[%d], y[%d]", x.length, y.length));
       }

       for (int i = 0; i < x.length; i++) {
           y[i] -= x[i];
       }
   }
   
   public static void main(String args[]) {
      MultivariateGaussianDistribution mvgd = new MultivariateGaussianDistribution(new double[] {20, 30}, MatrixUtils.createRealMatrix(new double[][]{{16., 21.6}, {21.6, 36.}}));
      System.out.println(mvgd.cdf(new double[] {10, 20}) - mvgd.cdf(new double[] {9, 19}));
   }
}
