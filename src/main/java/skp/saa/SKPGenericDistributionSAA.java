package skp.saa;

import java.util.Arrays;
import java.util.stream.IntStream;

import ilog.concert.IloException;
import skp.instance.SKPGenericDistribution;
import skp.milp.SKPScenarioBased;
import skp.milp.instance.SKPScenarioBasedSolvedInstance;
import skp.saa.instance.SKPGenericDistributionSAASolvedInstance;
import skp.utilities.gson.GSONUtility;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.rng.MRG32k3aL;
import umontreal.ssj.rng.RandomStream;

public class SKPGenericDistributionSAA {
   
   static final long[] seed = {12345, 24513, 24531, 42531, 35124, 32451};
   RandomStream randGenerator;
   
   SKPGenericDistribution instance;
   
   double optGap1;
   double optGap2;
   
   public SKPGenericDistributionSAA(SKPGenericDistribution instance) {
      MRG32k3aL rnd = new MRG32k3aL();
      rnd.setSeed(seed);
      this.randGenerator = rnd;
      this.randGenerator.resetStartStream();
      
      this.instance = instance;
   }
   
   /**
    * SAA estimators:
    * \hat{g}_{N'}(\hat{x})                  -> SKPScenarioBasedSolvedInstance.simulatedSolutionValueMean
    * S^2_{N'}(\hat{x})                      -> SKPScenarioBasedSolvedInstance.simulatedSolutionValueVariance
    * \hat{x}^M_N                            -> SKPScenarioBasedSolvedInstance.milpSolutionValue
    * \bar{v}^M_N                            -> barvMN
    * \hat{g}_{N'}(\hat{x}) - \bar{v}^M_N    -> optGapEstimator1
    * S^2_M/M                                -> barvMNVariance
    * S^2_{N'}(\hat{x}) + S^2_M/M            -> optGapEstimator1Variance
    * 
    * Alternative
    * \hat{g}^m_N(\hat{x})                   -> computeSampleAverage
    * \bar{g}^m_N(\hat{x}) - \bar{v}^M_N     -> optGapEstimator2     
    * \bar{S}^2_{M}/M                        -> optGapEstimator2Variance        
    */
   
   private static long time_limitMs = 60*10*1000; //10 minutes
   private static double tolerance = 1e-2;
   
   public SKPGenericDistributionSAASolvedInstance solve(int Nsmall, int Nlarge, int M) {
      long startGlobal = System.currentTimeMillis();
      SKPScenarioBasedSolvedInstance[] SAAreplications = new SKPScenarioBasedSolvedInstance[M];
      double barvMN = 0;
      
      double[][][] scenarios = new double[M][][]; 
      
      int bestSoFar = 0;
      double bestObjSoFar = Double.MIN_VALUE;
      for(int m = 0; m < M; m++) {
         try {
            int numberOfScenarios = Nsmall;
            SKPScenarioBased sskp = new SKPScenarioBased(instance, numberOfScenarios, this.randGenerator);
            scenarios[m] = sskp.scenarios;
            
            int simulationRuns = Nlarge;
            //System.out.println(GSONUtility.<SKPScenarioBasedSolvedInstance>printInstanceAsJSON(sskp.solve(simulationRuns)));
            SAAreplications[m] = sskp.solve(simulationRuns);
            barvMN += SAAreplications[m].milpSolutionValue;
            if(bestObjSoFar < SAAreplications[m].simulatedSolutionValueMean) {
               bestObjSoFar = SAAreplications[m].simulatedSolutionValueMean;
               bestSoFar = m;
            }
         } catch (IloException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
         }
         
         //Early stopping
         if(m > 32) { // 32 periods warm up for CLT to apply
            
            computeOptimalityGaps(SAAreplications, barvMN, scenarios, bestSoFar, bestObjSoFar, m + 1);
            
            
            
            if(this.optGap1 < tolerance || 
                  this.optGap2 < tolerance ||
                  System.currentTimeMillis() - startGlobal > time_limitMs) {
               double endGlobal = System.currentTimeMillis();
               double solutionTimeMs = endGlobal - startGlobal;
               
               return new SKPGenericDistributionSAASolvedInstance(this.instance, SAAreplications[bestSoFar].optimalKnapsack, SAAreplications[bestSoFar].simulatedSolutionValueMean, solutionTimeMs, this.optGap1, this.optGap2, Nsmall, Nlarge, m);
            }
         }
      }
      computeOptimalityGaps(SAAreplications, barvMN, scenarios, bestSoFar, bestObjSoFar, M);
      
      double endGlobal = System.currentTimeMillis();
      double solutionTimeMs = endGlobal - startGlobal;
      
      return new SKPGenericDistributionSAASolvedInstance(this.instance, SAAreplications[bestSoFar].optimalKnapsack, SAAreplications[bestSoFar].simulatedSolutionValueMean, solutionTimeMs, this.optGap1, this.optGap2, Nsmall, Nlarge, M);
   }

   private void computeOptimalityGaps(SKPScenarioBasedSolvedInstance[] SAAreplications, double barvMN,
         double[][][] scenarios, int bestSoFar, double bestObjSoFar, int W) {
      final double final_barvMN = barvMN/W;
      final int final_bestSoFar = bestSoFar;
      
      // Estimators
      double optGapEstimator1 = bestObjSoFar - final_barvMN;
      double barvMNVariance = 1.0/(W*(W-1))*Arrays.stream(SAAreplications).limit(W).mapToDouble(r -> Math.pow(r.milpSolutionValue - final_barvMN, 2)).sum();
      double optGapEstimator1Variance = SAAreplications[bestSoFar].simulatedSolutionValueVariance + barvMNVariance;
      
      double bargMN = IntStream.iterate(0, i -> i + 1).limit(W).mapToDouble(r -> computeSampleAverage(SAAreplications[final_bestSoFar].optimalKnapsack, scenarios[r])).sum()/W;
      double optGapEstimator2 = bargMN - final_barvMN;
      double optGapEstimator2Variance = 1.0/(W*(W-1))*IntStream.iterate(0, i -> i + 1).limit(W).mapToDouble(r -> Math.pow((computeSampleAverage(SAAreplications[final_bestSoFar].optimalKnapsack, scenarios[r])-SAAreplications[r].milpSolutionValue)-(bargMN - final_barvMN),2)).sum();
      
      double z = NormalDist.inverseF01(0.95);
      this.optGap1 = - (optGapEstimator1 - z * Math.sqrt(optGapEstimator1Variance)); // We are in a max setting!
      this.optGap2 = - (optGapEstimator2 - z * Math.sqrt(optGapEstimator2Variance)); // We are in a max setting!
   }
   
   public double computeSampleAverage(int[] knapsack, double[][] sampleMatrix) {
      double knapsackValue = 0;
      for(int i = 0; i < knapsack.length; i++) {
         if(knapsack[i] == 1) knapsackValue += this.instance.getExpectedValues()[i]; 
      }
      knapsackValue *= sampleMatrix.length;
      for(int s = 0; s < sampleMatrix.length; s++) {
         double weight = 0;
         for(int i = 0; i < knapsack.length; i++) {
            if(knapsack[i] == 1) {
               weight += sampleMatrix[s][i];
            }
         }
         knapsackValue -= instance.getShortageCost()*Math.max(0, weight - instance.getCapacity());
         
      }
      knapsackValue /= sampleMatrix.length;
      return knapsackValue;
   }
   
   public static void main(String args[]) {

      SKPGenericDistribution instance = SKPGenericDistribution.getTestInstanceLarge();
      
      SKPGenericDistributionSAA skpSAA = new SKPGenericDistributionSAA(instance);
      
      // SAA parameters
      int Nsmall = 100;    // N
      int Nlarge = 100000; // N'
      int M = 100;         // M
      
      System.out.println(GSONUtility.<SKPGenericDistributionSAASolvedInstance>printInstanceAsJSON(skpSAA.solve(Nsmall, Nlarge, M)));
   }
}
