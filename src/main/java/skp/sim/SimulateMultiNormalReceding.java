/**
 * To run from Mac OS
 * 
 * -Djava.library.path=/Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
 * 
 * Environment variable
 * 
 * DYLD_LIBRARY_PATH /Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
 * 
 * @author gwren
 *
 */

package skp.sim;

import java.util.Arrays;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import ilog.concert.IloException;
import skp.KP;
import skp.SKPMultiNormal;
import skp.milp.KPMILP;
import skp.milp.SKPMultiNormalMILP;
import umontreal.ssj.probdistmulti.MultiNormalDist;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalCholeskyGen;
import umontreal.ssj.randvarmulti.MultinormalGen;
import umontreal.ssj.rng.MRG32k3aL;

/**
 * TO DO: Remember to amend representation of item value (per unit and proportional to item weight)
 * 
 * @author gwren
 *
 */

public class SimulateMultiNormalReceding {
   SKPMultiNormal instance;
   private MRG32k3aL randGenerator;
   
   public SimulateMultiNormalReceding(SKPMultiNormal instance, long[] seed) {
      this.randGenerator = new MRG32k3aL();
      this.randGenerator.setSeed(seed);
      this.instance = instance;
   }
   
   private static RealMatrix matrixInverse(RealMatrix m) {
      RealMatrix pInverse = new LUDecomposition(m).getSolver().getInverse();
      return pInverse;
   }
   
   private static RealMatrix createReducedMatrix(RealMatrix m) {
      double[][] matrix = new double[m.getRowDimension()-1][m.getColumnDimension()-1];
      for(int i = 1; i < m.getRowDimension(); i++) {
         for(int j = 1; j < m.getColumnDimension(); j++) {
            matrix[i-1][j-1]=m.getEntry(i, j);
         }
      }
      RealMatrix n = MatrixUtils.createRealMatrix(matrix);
      return n;  
   }
   
   private static RealMatrix computeConditionalExpectedDemand(double [] expDemand, double [] realizationDemand, RealMatrix m) {
      RealMatrix d1 =null;
      RealMatrix d2 =null;
      RealMatrix sigma21 =null;
      RealMatrix sigma11Inv =null;
      RealMatrix zeta1 =null;
      
      
      double[][] ed1 = new double[1][];
      ed1[0] = expDemand;
      RealMatrix ed1Matrix = MatrixUtils.createRealMatrix(ed1); 
      d1=ed1Matrix.getSubMatrix(0, 0, 0, 0).transpose(); //d1
      
   
      double[][] ed2 = new double[1][];
      ed2[0] = expDemand;
      RealMatrix ed2Matrix = MatrixUtils.createRealMatrix(ed2); 
      d2=ed2Matrix.getSubMatrix(0, 0, 1, expDemand.length - 1).transpose(); //d2
      
      sigma21 = m.getSubMatrix(1, m.getRowDimension() - 1, 0, 0);
      
      sigma11Inv = matrixInverse(m.getSubMatrix(0, 0, 0, 0));
      
      double[][] realisedDemand = new double[1][];
      realisedDemand[0] = realizationDemand;
      RealMatrix realisedDemandMatrix = MatrixUtils.createRealMatrix(realisedDemand); 
      zeta1 = realisedDemandMatrix.getSubMatrix(0, 0, 0, 0).transpose();
      
      
      RealMatrix results = d2.add(sigma21.multiply(sigma11Inv.multiply(zeta1.subtract(d1))));      
      
        return results.transpose();    
   }
   
   private int simulateOneItem(int t, double[] realizations, double remainingCapacity, int partitions) {
      double[] previousRealizations = new double[t];
      System.arraycopy(realizations, 0, previousRealizations, 0, t);
      
      double[] expDemand = instance.getWeights().getMean();
      
      double[][] covariance = instance.getWeights().getCovariance();
      
      int Nbperiods = expDemand.length-t;
      
      double[] shortexpValues = new double[Nbperiods];
      System.arraycopy(instance.getExpectedValues(), t, shortexpValues, 0, Nbperiods);
      double[] shortexpDemand = expDemand;
      RealMatrix conditionalCovarianceMatrix = MatrixUtils.createRealMatrix(covariance);
      RealMatrix M = conditionalCovarianceMatrix;
      double[] realizationsBuffer = realizations;
      for(int k = 0; k < t; k++) {
         RealMatrix inverseM =  matrixInverse(conditionalCovarianceMatrix);
         RealMatrix reducedInverseM = createReducedMatrix(inverseM);
         conditionalCovarianceMatrix =  matrixInverse(reducedInverseM);
         shortexpDemand = computeConditionalExpectedDemand(shortexpDemand, realizationsBuffer, M).getRow(0);
         double[][] realizationMatrix = new double[1][];
         realizationMatrix[0] = realizationsBuffer;
         RealMatrix reducedRealizations = MatrixUtils.createRealMatrix(realizationMatrix);
         realizationsBuffer =  reducedRealizations.getSubMatrix(0, 0, 1, realizationsBuffer.length - 1).getRow(0);
         M = conditionalCovarianceMatrix;
      }
      double [][] shortCovariance = conditionalCovarianceMatrix.getData();
      
      MultiNormalDist weights = new MultiNormalDist(shortexpDemand, shortCovariance);
      
      SKPMultiNormal reducedInstance = new SKPMultiNormal(shortexpValues, weights, remainingCapacity, instance.getShortageCost());
      
      SKPMultiNormalMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new SKPMultiNormalMILP(reducedInstance, partitions);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      
      return knapsack[0];
   }
   
   private double simulateOneRun(double[] realizations, int partitions) {
      double knapsackValue = 0;
      double remainingCapacity = instance.getCapacity();
      for(int i = 0; i < realizations.length; i++) {
         if(simulateOneItem(i, realizations, remainingCapacity, partitions) == 1) {
            remainingCapacity -= realizations[i];
            knapsackValue += this.instance.getExpectedValues()[i];
         }
      }
      knapsackValue -= Math.max(-remainingCapacity*instance.getShortageCost(), 0);
      return knapsackValue;
   }
   
   public double simulate(int nbSamples, int partitions) {
      double[][] sampleMatrix = sampleWeights(nbSamples);
      double knapsackValue = Arrays.stream(sampleMatrix).parallel().mapToDouble(r -> simulateOneRun(r, partitions)).sum()/nbSamples;
      return knapsackValue;
   }
   
   private double simulateOneRunEVPI(double[] realizations) {
      KP instanceEVPI = new KP(instance.getExpectedValues(), realizations, instance.getCapacity(), instance.getShortageCost());
      KPMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new KPMILP(instanceEVPI);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      return milp.getSolutionValue();
   }
   
   public double simulateEVPI(int nbSamples) {
      double[][] sampleMatrix = sampleWeights(nbSamples);
      double knapsackValue = Arrays.stream(sampleMatrix).parallel().mapToDouble(r -> simulateOneRunEVPI(r)).sum()/nbSamples;
      return knapsackValue;
   }
   
   public double computeEVP() {
      KP instanceEVP = new KP(instance.getExpectedValues(), instance.getWeights().getMean(), instance.getCapacity(), instance.getShortageCost());
      KPMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new KPMILP(instanceEVP);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      return milp.getSolutionValue();
   }
   
   private double simulateOneRunEVwPI(double[] realizations, int selectFirstObject) {
      KP instanceEVPI = new KP(instance.getExpectedValues(), realizations, instance.getCapacity(), instance.getShortageCost());
      KPMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new KPMILP(instanceEVPI, selectFirstObject);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      return milp.getSolutionValue();
   }
   
   public double simulateEVwPI(int nbSamples) {
      double[][] sampleMatrix = sampleWeights(nbSamples);
      double knapsackValueX1Eq0 = Arrays.stream(sampleMatrix).parallel().mapToDouble(r -> simulateOneRunEVwPI(r, 0)).sum()/nbSamples;
      double knapsackValueX1Eq1 = Arrays.stream(sampleMatrix).parallel().mapToDouble(r -> simulateOneRunEVwPI(r, 1)).sum()/nbSamples;
      return Math.max(knapsackValueX1Eq0, knapsackValueX1Eq1);
   }

   
   private double[][] sampleWeights(int nbSamples) {
      double[] mu = this.instance.getWeights().getMean();
      double[][] sigma = this.instance.getWeights().getCovariance();
      this.randGenerator.resetStartStream();
      NormalGen standardNormal = new NormalGen(this.randGenerator, 0, 1);
      MultinormalGen gen = new MultinormalCholeskyGen(standardNormal, mu, sigma);
      double[][] points = new double[nbSamples][this.instance.getItems()];
      for(int i = 0; i < nbSamples; i++)
         gen.nextPoint(points[i]);
      return points;
   }
   
   private static double[][] calculateCovariance(double [] means, double cv, double rho){
      double [] stdDemand =new double [means.length];
      for (int i = 0; i < means.length; i ++) {
         stdDemand[i] = cv*means[i];
      }
      
      double [][] covariance = new double [means.length][means.length];
      
      for (int row=0; row<covariance.length;row++) {
         for (int col=0; col<covariance[row].length;col++) {
            if (row==col) {
               covariance[row][col]=stdDemand[row]*stdDemand[col];
            } else if (col==row+1 | col==row-1) {
               covariance[row][col]=stdDemand[row]*stdDemand[col]*rho;
            } else  {
               covariance[row][col]=0;
            }
         }
      }
      return covariance;
   }
   
   public static void main(String args[]) {
      long[] seed = {1,2,3,4,5,6};
      
      double[] expectedValues = {111,111,21,117,123,34,3,121,112,12};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      double rho = 0.5;
      double[][] varianceCovarianceWeights = calculateCovariance(expectedWeights, cv, rho);
      
      int capacity = 100;
      int shortageCost = 100;
      
      MultiNormalDist weights = new MultiNormalDist(expectedWeights, varianceCovarianceWeights);
      SKPMultiNormal instance = new SKPMultiNormal(expectedValues, weights, capacity, shortageCost);
      
      int partitions = 10;
      int nbSamples = 20;
      
      SimulateMultiNormalReceding sim = new SimulateMultiNormalReceding(instance, seed);
      
      double simSolutionValue = sim.simulate(nbSamples, partitions);
      double simEVPI = sim.simulateEVPI(nbSamples);
      double EVP = sim.computeEVP();
      double EVwP = sim.simulateEVwPI(nbSamples);
      
      System.out.println("Simulation: "+simSolutionValue);
      System.out.println("EVPI: "+simEVPI);
      System.out.println("EVP: "+EVP);
      System.out.println("EVwP on items 2-n: "+EVwP);
   }
}


