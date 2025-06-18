package skp.folf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.PoissonDist;

public class PiecewisePoissonFirstOrderLossFunction {

   /**
    * Poisson knapsack
    */
   
   /**
    *  Memoization
    */
   static Hashtable<String,double[]> poisson_ce = new Hashtable<String,double[]>();
   static Hashtable<String,Double> poisson_ae = new Hashtable<String,Double>();
   
   public static double[][] poissonKnapsackPiecewiseFOLFConditionalExpectations(int capacity, double[] probabilityMasses, int nbSamples) {
      ArrayList<double[]> conditionalExpectations = new ArrayList<double[]>();
      ArrayList<Double> approximationErrors = new ArrayList<Double>();
      
      int demand = 0;
      while(true) {
         String hash = demand + "_" + Arrays.toString(probabilityMasses) + "_" + nbSamples;
         Distribution[] distributions = new Distribution[1];
         distributions[0] = new PoissonDist(demand);
         PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
         double[] ce = null;
         if(poisson_ce.containsKey(hash)) {
            ce = poisson_ce.get(hash);
         } else {
            ce = pwfolf.getConditionalExpectations(probabilityMasses, nbSamples);
            poisson_ce.put(hash, ce);
         }
         conditionalExpectations.add(ce);
         double ae;
         if(poisson_ae.containsKey(hash)) {
            ae = poisson_ae.get(hash);
         } else {
            ae = pwfolf.getMaxApproximationError(probabilityMasses, nbSamples);
            poisson_ae.put(hash, ae);
         }
         approximationErrors.add(ae);
         if(ce[0] > capacity) break;
         else demand++;
      }
      return conditionalExpectations.toArray(new double[conditionalExpectations.size()][]);
   }
   public static double[] poissonKnapsackPiecewiseFOLFApproximationErrors(int capacity, double[] probabilityMasses, int nbSamples) {
      ArrayList<double[]> conditionalExpectations = new ArrayList<double[]>();
      ArrayList<Double> approximationErrors = new ArrayList<Double>();
      
      int demand = 0;
      while(true) {
         String hash = demand + "_" + Arrays.toString(probabilityMasses) + "_" + nbSamples;
         Distribution[] distributions = new Distribution[1];
         distributions[0] = new PoissonDist(demand);
         PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
         double[] ce = null;
         if(poisson_ce.containsKey(hash)) {
            ce = poisson_ce.get(hash);
         } else {
            ce = pwfolf.getConditionalExpectations(probabilityMasses, nbSamples);
            poisson_ce.put(hash, ce);
         }
         conditionalExpectations.add(ce);
         double ae;
         if(poisson_ae.containsKey(hash)) {
            ae = poisson_ae.get(hash);
         } else {
            ae = pwfolf.getMaxApproximationError(probabilityMasses, nbSamples);
            poisson_ae.put(hash, ae);
         }
         approximationErrors.add(ae);
         if(ce[0] > capacity) break;
         else demand++;
      }
      return PiecewiseFirstOrderLossFunction.toPrimitive(approximationErrors.toArray(new Double[conditionalExpectations.size()]));
   }
   public static int poissonKnapsackPiecewiseFOLFMaxWeight(int capacity, double[] probabilityMasses, int nbSamples) {
      ArrayList<double[]> conditionalExpectations = new ArrayList<double[]>();
      ArrayList<Double> approximationErrors = new ArrayList<Double>();
      
      int demand = 0;
      while(true) {
         String hash = demand + "_" + Arrays.toString(probabilityMasses) + "_" + nbSamples;
         Distribution[] distributions = new Distribution[1];
         distributions[0] = new PoissonDist(demand);
         PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
         double[] ce = null;
         if(poisson_ce.containsKey(hash)) {
            ce = poisson_ce.get(hash);
         } else {
            ce = pwfolf.getConditionalExpectations(probabilityMasses, nbSamples);
            poisson_ce.put(hash, ce);
         }
         conditionalExpectations.add(ce);
         double ae;
         if(poisson_ae.containsKey(hash)) {
            ae = poisson_ae.get(hash);
         } else {            
            ae = pwfolf.getMaxApproximationError(probabilityMasses, nbSamples);
            poisson_ae.put(hash, ae);
         }
         approximationErrors.add(ae);
         if(ce[0] > capacity) break;
         else demand++;
      }
      return demand;
   }
   @SuppressWarnings("unused")
   private static void generatePoissonPiecewiseTables() {
      int capacity = 100;
      int partitions = 5;
      int nbSamples = 10000;
      double[] probabilityMasses = new double[partitions];
      Arrays.fill(probabilityMasses, 1.0/partitions);
      System.out.println("prob = " + Arrays.toString(probabilityMasses) + ";");
      System.out.println("means = " + Arrays.deepToString(poissonKnapsackPiecewiseFOLFConditionalExpectations(capacity, probabilityMasses, nbSamples)) + ";");
      System.out.println("error = " + Arrays.toString(poissonKnapsackPiecewiseFOLFApproximationErrors(capacity, probabilityMasses, nbSamples)) + ";");
      System.out.println("maxWeight = " + poissonKnapsackPiecewiseFOLFMaxWeight(capacity, probabilityMasses, nbSamples) + ";");
   }
   @SuppressWarnings("unused")
   private static void printLinearizationParameters() {
      Distribution[] distributions = new Distribution[3];
      double lambda[] = {20,5,50};
      distributions[0] = new PoissonDist(lambda[0]);
      distributions[1] = new PoissonDist(lambda[1]);
      distributions[2] = new PoissonDist(lambda[2]);
      PiecewiseFirstOrderLossFunction pwfolf = new PiecewiseFirstOrderLossFunction(distributions);
      
      int partitions = 5;
      int nbSamples = 1000;
      double[] probabilityMasses = new double[partitions]; 
      Arrays.fill(probabilityMasses, 1.0/partitions);
      String out = "float prob[partitions] = ";
      out += Arrays.toString(probabilityMasses) + ";\n";
      out += "float means[partitions] = ";
      out += Arrays.toString(pwfolf.getConditionalExpectations(probabilityMasses, nbSamples)) + ";\n";
      out += "float error = "+pwfolf.getMaxApproximationError(probabilityMasses, nbSamples) + ";\n";
      System.out.print(out);
   }
   
   /**
    * Main
    */
   
   public static void main(String[] args){
      //printLinearizationParameters();
      //generatePoissonPiecewiseTables();
   }

}
