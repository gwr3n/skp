package skp.sdp;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import com.sun.management.OperatingSystemMXBean;

import gnu.trove.map.hash.THashMap;

import skp.instance.SKP;
import skp.instance.SKPMultinormal;
import skp.sdp.instance.DSKPMultinormalSolvedInstance;
import skp.utilities.gson.GSONUtility;

import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdistmulti.MultiNormalDist;

public class DSKPMultinormal{
   
   int[] supportLB;
   int[] supportUB;
   
   OperatingSystemMXBean osMBean;
   long nanoBefore;
   long nanoAfter;
   
   Function<State, double[]> actionGenerator;
   
   @FunctionalInterface
   interface StateTransitionFunction <S, M, A> { 
      public S apply (S s, M m, A a);
   }
   
   public StateTransitionFunction<State, Double, Double> stateTransition;
   
   @FunctionalInterface
   interface ImmediateValueFunction <S, A, R, V> { 
      public V apply (S s, A a, R r);
   }
   
   public ImmediateValueFunction<State, Double, Double, Double> immediateValueFunction;
   
   public void initialiseFunctionalInterfaces(SKP instance) {
      /**
       * This function returns the set of actions associated with a given state
       */
      this.actionGenerator = state ->{
         return DoubleStream.iterate(0, num -> num + 1)
                            .limit(2)
                            .toArray();
      };
      
      /**
       * State transition function; given a state, an action and a random outcome, the function
       * returns the future state
       */
      this.stateTransition = (state, previousItemWeight, realisedWeight) -> 
         this.new State(state.item + 1, previousItemWeight, state.remainingCapacity - realisedWeight);
      
      /**
       * Immediate value function for a given state
       */
      this.immediateValueFunction = (state, action, realisedWeight) -> {
            double value = action*instance.getExpectedValuesPerUnit()[state.item]*realisedWeight;
            double cost = (state.item == instance.getItems() - 1 ? instance.getShortageCost() : 0)*Math.max(action*realisedWeight - state.remainingCapacity, 0);
            return value - cost;
         };
   }
   
   Map<State, Double> cacheActions = Collections.synchronizedMap(new THashMap<>());
   Map<State, Double> cacheValueFunction = Collections.synchronizedMap(new THashMap<>());
   
   class State{
      int item;  
      double previousItemWeight;
      double remainingCapacity;            
      
      public State(int item, 
                   double previousItemWeight,
                   double remainingCapacity){
         this.item = item;
         this.previousItemWeight = previousItemWeight;
         this.remainingCapacity = remainingCapacity;
      }
      
      public double[] getFeasibleActions(){
         return actionGenerator.apply(this);
      }
      
      @Override
      public int hashCode(){
         String hash = "";
         hash = (hash + this.item) + "_" + this.previousItemWeight + "_" + this.remainingCapacity;
         return hash.hashCode();
      }
      
      @Override
      public boolean equals(Object o){
         if(o instanceof State)
            return  ((State) o).item == this.item &&
                    ((State) o).previousItemWeight == this.previousItemWeight &&
                    ((State) o).remainingCapacity == this.remainingCapacity;
         else
            return false;
      }

      @Override
      public String toString(){
         return this.item + " " + this.previousItemWeight + " " + this.remainingCapacity;
      }
   }
   
   SKPMultinormal instance;
   static double d = 0.5;
   
   /**
    * Computes the optimal solution to the SKPMultinormal <code>instance</code> using Stochastic Dynamic Programming.<br>
    * <br>
    * Note that the weight covariance matrix of <code>instance</code> must be generated according to the special 
    * structure $\rho^{|j-i|}\sigma_i\sigma_j$ discussed in [1], which ensures $P(d_t=x|d_{t-1}=y) = P(d_t=x|d_{t-1}=y,d_{t-2}=z,...)$.<br>
    * <br>
    * [1]M. Xiang, R. Rossi, B. Martin-Barragan, S. A. Tarim, "<a href="https://doi.org/10.1016/j.ejor.2022.04.011">
    * A mathematical programming-based solution method for the nonstationary inventory problem under correlated demand</a>,
    * " European Journal of Operational Research, Elsevier, Vol. 304(2): 515–524, 2023
    * 
    */
   public DSKPMultinormal(SKPMultinormal instance, double truncationQuantile) {
      this.instance = instance;
      supportLB = IntStream.iterate(0, i -> i + 1)
                           .limit(instance.getWeights().getMu().length)
                           .map(i -> (int)Math.round((new NormalDist(instance.getWeights().getMu(i),Math.sqrt(instance.getWeights().getSigma()[i][i]))).inverseF(1-truncationQuantile)))
                           .toArray();
      supportUB = IntStream.iterate(0, i -> i + 1)
                           .limit(instance.getWeights().getMu().length)
                           .map(i -> (int)Math.round((new NormalDist(instance.getWeights().getMu(i),Math.sqrt(instance.getWeights().getSigma()[i][i]))).inverseF(truncationQuantile)))
                           .toArray();
      
      initialiseFunctionalInterfaces(instance);
   }
   
   public DSKPMultinormalSolvedInstance solve() {
      /**
       * Initial problem conditions
       */
      int initialItem = 0;
      State initialState = this.new State(initialItem, Double.NaN, instance.getCapacity());
      
      try {
         osMBean = ManagementFactory.newPlatformMXBeanProxy(
               ManagementFactory.getPlatformMBeanServer(), ManagementFactory.OPERATING_SYSTEM_MXBEAN_NAME, OperatingSystemMXBean.class);
         nanoBefore = System.nanoTime();
      } catch (IOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
      
      double solutionValue = this.f(initialState);
      this.nanoAfter = System.nanoTime();
      long solutionTimeMs = (int) Math.ceil(((this.nanoAfter-this.nanoBefore)*Math.pow(10, -6)));
      long statesExplored = this.cacheValueFunction.size();
      
      return new DSKPMultinormalSolvedInstance(instance, solutionValue, solutionTimeMs, statesExplored);
   }
   
   double f(State state){
      //if(state.item == 1) System.out.println(state);
      return cacheValueFunction.computeIfAbsent(state, s -> {
         double val= Arrays.stream(s.getFeasibleActions())
                           .map(x ->          DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                          .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                          .map(w -> getMarginalCDFDifference(instance.getWeights(), s.item, s.previousItemWeight, w)*immediateValueFunction.apply(s, x, w)+ 
                                                                    ((s.item < instance.getItems() - 1) ? getMarginalCDFDifference(instance.getWeights(), s.item, s.previousItemWeight, w)*f(stateTransition.apply(s, w, x == 0 ? 0 : w)) : 0))
                                                          .sum())
                           .max()
                           .getAsDouble();
         double bestAction = Arrays.stream(s.getFeasibleActions())
                                .filter(x -> (DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                          .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                          .map(w -> getMarginalCDFDifference(instance.getWeights(), s.item, s.previousItemWeight, w)*immediateValueFunction.apply(s, x, w)+ 
                                                                    ((s.item < instance.getItems() - 1) ? getMarginalCDFDifference(instance.getWeights(), s.item, s.previousItemWeight, w)*f(stateTransition.apply(s, w, x == 0 ? 0 : w)) : 0))
                                                          .sum()) == val)
                                .findAny()
                                .getAsDouble();
         cacheActions.putIfAbsent(s, bestAction);
         return val;
      });
   }
   
   public static void main(String args[]) {

      SKPMultinormal instance = SKPMultinormal.getTestInstanceSpecialStructure();
      
      double truncationQuantile = 0.9999;
      
      DSKPMultinormal dskp = new DSKPMultinormal(instance, truncationQuantile);
      
      System.out.println(GSONUtility.<DSKPMultinormalSolvedInstance>printInstanceAsJSON(dskp.solve()));
   }
   
   public static double getMarginalCDFDifference(MultiNormalDist dist, int coordinate, double realisedDemand, double value) {
      if(coordinate == 0) {
         NormalDist normal = new NormalDist(dist.getMu(0), Math.sqrt(dist.getSigma()[0][0]));
         return normal.cdf(value + d) - normal.cdf(value - d);
      } else {
         double[][] conditionalCovarianceMatrix =  MatrixAlgebra.computeConditionalCovarianceMatrix(dist.getCovariance(), coordinate);
         double[] realizationsBuffer = new double[coordinate];
         realizationsBuffer[coordinate-1] = realisedDemand;
         double[] shortExpWeight = MatrixAlgebra.computeConditionalExpectedDemand(dist.getMu(), realizationsBuffer, MatrixUtils.createRealMatrix(dist.getCovariance()));
         //MultivariateGaussianDistribution mvgd = new MultivariateGaussianDistribution(shortExpWeight, MatrixUtils.createRealMatrix(conditionalCovarianceMatrix));
         //return mvgd.cdf(new double[] {value + d}) - mvgd.cdf(new double[] {value - d});
         NormalDist normal = new NormalDist(shortExpWeight[0], Math.sqrt(conditionalCovarianceMatrix[0][0]));
         return normal.cdf(value + d) - normal.cdf(value - d);
      }
   }
   
private static class MatrixAlgebra {
      
      private static RealMatrix matrixInverse(RealMatrix m) {
         RealMatrix pInverse = new LUDecomposition(m).getSolver().getInverse();
         return pInverse;
      }
      
      private static RealMatrix createReducedMatrix(RealMatrix m, int period) {
         double[][] matrix = new double[m.getRowDimension() - period][m.getColumnDimension( )- period];
         for(int i = period; i < m.getRowDimension(); i++) {
            for(int j = period; j < m.getColumnDimension(); j++) {
               matrix[i-period][j-period]=m.getEntry(i, j);
            }
         }
         RealMatrix n = MatrixUtils.createRealMatrix(matrix);
         return n;  
      }
      
      public static double[][] computeConditionalCovarianceMatrix(double[][] covarianceMatrix, int period) {
         RealMatrix conditionalCovarianceMatrix = MatrixUtils.createRealMatrix(covarianceMatrix);
         RealMatrix inverseM =  MatrixAlgebra.matrixInverse(conditionalCovarianceMatrix);
         RealMatrix reducedInverseM = MatrixAlgebra.createReducedMatrix(inverseM, period);
         return MatrixAlgebra.matrixInverse(reducedInverseM).getData();
      }
      
      public static double[] computeConditionalExpectedDemand(double [] expDemand, double [] realizationDemand, RealMatrix m) {
         RealMatrix d1 =null;
         RealMatrix d2 =null;
         RealMatrix sigma21 =null;
         RealMatrix sigma11Inv =null;
         RealMatrix zeta1 =null;
         
         int period = realizationDemand.length;
         
         double[][] ed1 = new double[1][];
         ed1[0] = expDemand;
         RealMatrix ed1Matrix = MatrixUtils.createRealMatrix(ed1); 
         d1=ed1Matrix.getSubMatrix(0, 0, 0, period - 1).transpose(); //d1
         
         double[][] ed2 = new double[1][];
         ed2[0] = expDemand;
         RealMatrix ed2Matrix = MatrixUtils.createRealMatrix(ed2); 
         d2=ed2Matrix.getSubMatrix(0, 0, period, expDemand.length - 1).transpose(); //d2
         
         sigma21 = m.getSubMatrix(period, m.getRowDimension() - 1, 0, period - 1);
         
         sigma11Inv = matrixInverse(m.getSubMatrix(0, period - 1, 0, period - 1));
         
         double[][] realisedDemand = new double[1][];
         realisedDemand[0] = realizationDemand;
         RealMatrix realisedDemandMatrix = MatrixUtils.createRealMatrix(realisedDemand); 
         zeta1 = realisedDemandMatrix.transpose();
         
         RealMatrix results = d2.add(sigma21.multiply(sigma11Inv.multiply(zeta1.subtract(d1))));      
         
         return results.transpose().getRow(0);    
      }
   }
}


