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
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import ilog.concert.IloException;

import skp.folf.PiecewiseStandardNormalFirstOrderLossFunction;
import skp.instance.KP;
import skp.instance.SKPMultinormal;
import skp.milp.KPMILP;
import skp.milp.SKPMultinormalMILP;
import skp.sim.instance.SKPMultinormalRecedingSolvedInstance;

import umontreal.ssj.probdistmulti.MultiNormalDist;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalGen;
import umontreal.ssj.randvarmulti.MultinormalPCAGen;

public class SimulateMultinormalReceding  extends Simulate {
   
   SKPMultinormal instance;
   int partitions;
   
   public SimulateMultinormalReceding(SKPMultinormal instance, int partitions) {
      this.instance = instance;
      this.partitions = partitions;
   }
   
   public SKPMultinormalRecedingSolvedInstance solve(int simulationRuns) {
      
      SimulateMultinormalReceding sim = new SimulateMultinormalReceding(instance, partitions);
      
      double[] realisations = sim.simulate(simulationRuns, partitions);
      Mean m = new Mean();
      double simSolutionMean = m.evaluate(realisations);
      StandardDeviation std = new StandardDeviation();
      double simSolutionStd = std.evaluate(realisations);
      double simEVwPI = sim.simulateEVwPI(simulationRuns);
      double EVP = sim.computeEVP();
      double EVwPI_obj_2_n = sim.simulateEVwPI_obj_2_n(simulationRuns);
      
      System.out.println("Simulation (mean): "+simSolutionMean);
      System.out.println("Simulation (std): "+simSolutionStd);
      System.out.println("EVwPI: "+simEVwPI);
      System.out.println("EVP: "+EVP);
      System.out.println("EVwPI on items 2-n: "+EVwPI_obj_2_n);
      
      SKPMultinormalRecedingSolvedInstance solvedInstance = new SKPMultinormalRecedingSolvedInstance(
            this.instance,
            simSolutionMean,
            simSolutionStd,
            simulationRuns,
            EVP,
            simEVwPI,
            EVwPI_obj_2_n,
            partitions,
            PiecewiseStandardNormalFirstOrderLossFunction.getLinearizationSamples()
            );
      
      return solvedInstance;
   }
   
   private int simulateOneItem(int t, double[] realizations, double remainingCapacity, int partitions) {
      double[] previousRealizations = new double[t];
      System.arraycopy(realizations, 0, previousRealizations, 0, t);
      
      double[] expWeight = instance.getWeights().getMean();
      double[][] covariance = instance.getWeights().getCovariance();
      
      int Nbperiods = expWeight.length-t;
      
      double[] shortExpValues = new double[Nbperiods];
      System.arraycopy(instance.getExpectedValuesPerUnit(), t, shortExpValues, 0, Nbperiods);
      
      double[] shortExpWeight = expWeight;
      RealMatrix conditionalCovarianceMatrix = MatrixUtils.createRealMatrix(covariance);
      RealMatrix M = conditionalCovarianceMatrix;
      double[] realizationsBuffer = realizations;
      for(int k = 0; k < t; k++) {
         RealMatrix inverseM =  MatrixAlgebra.matrixInverse(conditionalCovarianceMatrix);
         RealMatrix reducedInverseM = MatrixAlgebra.createReducedMatrix(inverseM);
         conditionalCovarianceMatrix =  MatrixAlgebra.matrixInverse(reducedInverseM);
         shortExpWeight = MatrixAlgebra.computeConditionalExpectedDemand(shortExpWeight, realizationsBuffer, M).getRow(0);
         double[][] realizationMatrix = new double[1][];
         realizationMatrix[0] = realizationsBuffer;
         RealMatrix reducedRealizations = MatrixUtils.createRealMatrix(realizationMatrix);
         realizationsBuffer =  reducedRealizations.getSubMatrix(0, 0, 1, realizationsBuffer.length - 1).getRow(0);
         M = conditionalCovarianceMatrix;
      }
      double [][] shortCovariance = conditionalCovarianceMatrix.getData();
      
      MultiNormalDist weights = new MultiNormalDist(shortExpWeight, shortCovariance);
      SKPMultinormal reducedInstance = new SKPMultinormal(shortExpValues, weights, remainingCapacity, instance.getShortageCost());
      
      SKPMultinormalMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new SKPMultinormalMILP(reducedInstance, partitions);
         milp.solve();
         knapsack = milp.getOptimalKnapsack();
         //System.out.println("Knapsack: "+Arrays.toString(knapsack));
         System.out.print(".");
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
            knapsackValue += this.instance.getExpectedValuesPerUnit()[i]*realizations[i];
         }
      }
      knapsackValue -= Math.max(-remainingCapacity*instance.getShortageCost(), 0);
      return knapsackValue;
   }
   
   double[] simulate(int nbSamples, int partitions) {
      double[][] sampleMatrix = sampleWeights(nbSamples);
      double[] knapsackValues = Arrays.stream(sampleMatrix)
                                      .parallel()
                                      .mapToDouble(r -> simulateOneRun(r, partitions))
                                      .peek(r -> System.out.println("Simulation run completed: "+r))
                                      .toArray();
      return knapsackValues;
   }
   
   private double simulateOneRunEVwPI(double[] realizations) {
      KP instanceEVwPI = new KP(instance.getExpectedValuesPerUnit(), realizations, instance.getCapacity(), instance.getShortageCost());
      KPMILP milp = null;
      int[] knapsack = null;
      try {
         milp = new KPMILP(instanceEVwPI);
         knapsack = milp.getKnapsack();
         System.out.println("Knapsack: "+Arrays.toString(knapsack));
      } catch (IloException e) {
         e.printStackTrace();
         System.exit(-1);
      }
      return milp.getSolutionValue();
   }
   
   double simulateEVwPI(int nbSamples) {
      double[][] sampleMatrix = sampleWeights(nbSamples);
      double knapsackValue = Arrays.stream(sampleMatrix).parallel().mapToDouble(r -> simulateOneRunEVwPI(r)).sum()/nbSamples;
      return knapsackValue;
   }
   
   double computeEVP() {
      KP instanceEVP = new KP(instance.getExpectedValuesPerUnit(), instance.getWeights().getMean(), instance.getCapacity(), instance.getShortageCost());
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
      KP instanceEVPI = new KP(instance.getExpectedValuesPerUnit(), realizations, instance.getCapacity(), instance.getShortageCost());
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
   
   double simulateEVwPI_obj_2_n(int nbSamples) {
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
      //MultinormalGen gen = new MultinormalCholeskyGen(standardNormal, mu, sigma);
      MultinormalGen gen = new MultinormalPCAGen(standardNormal, mu, sigma);
      double[][] points = new double[nbSamples][this.instance.getItems()];
      for(int i = 0; i < nbSamples; i++)
         gen.nextPoint(points[i]);
      return points;
   }
   
   public static void main(String args[]) {

      SKPMultinormal instance = SKPMultinormal.getTestInstanceSpecialStructure();
      
      int partitions = 10;
      int simulationRuns = 100;
      
      SimulateMultinormalReceding sim = new SimulateMultinormalReceding(instance, partitions);
      
      double[] realisations = sim.simulate(simulationRuns, partitions);
      Mean m = new Mean();
      double simSolutionMean = m.evaluate(realisations);
      StandardDeviation std = new StandardDeviation();
      double simSolutionStd = std.evaluate(realisations);
      double simEVwPI = sim.simulateEVwPI(simulationRuns);
      double EVP = sim.computeEVP();
      double EVwPI_obj_2_n = sim.simulateEVwPI_obj_2_n(simulationRuns);
      
      System.out.println("Simulation (mean): "+simSolutionMean);
      System.out.println("Simulation (std): "+simSolutionStd);
      System.out.println("Simulation (CI): ["+(simSolutionMean-1.96*simSolutionStd/Math.sqrt(simulationRuns))+","+
                                              (simSolutionMean+1.96*simSolutionStd/Math.sqrt(simulationRuns))+"]");
      System.out.println("EVwPI: "+simEVwPI);
      System.out.println("EVP: "+EVP);
      System.out.println("EVwPI on items 2-n: "+EVwPI_obj_2_n);
   }
   
   private static class MatrixAlgebra {
      
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
   }
}


