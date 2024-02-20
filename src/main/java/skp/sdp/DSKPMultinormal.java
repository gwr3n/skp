package skp.sdp;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;

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
   interface StateTransitionFunction <S, A> { 
      public S apply (S s, A a);
   }
   
   public StateTransitionFunction<State, Double> stateTransition;
   
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
      this.stateTransition = (state, realisedWeight) -> 
         this.new State(state.item + 1, realisedWeight, state.remainingCapacity - realisedWeight);
      
      /**
       * Immediate value function for a given state
       */
      this.immediateValueFunction = (state, action, realisedWeight) -> {
            double value = action*instance.getExpectedValuesPerUnit()[state.item]*realisedWeight;
            double cost = (state.item == instance.getItems() - 1 ? instance.getShortageCost() : 0)*Math.max(realisedWeight - state.remainingCapacity, 0);
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
      if(state.item == 1) System.out.println(state);
      return cacheValueFunction.computeIfAbsent(state, s -> {
         double val= Arrays.stream(s.getFeasibleActions())
                           .map(x -> (x == 0) ? immediateValueFunction.apply(s, 0.0, 0.0) + ((s.item < instance.getItems() - 1) ? f(stateTransition.apply(s, 0.0)) : 0) :
                                                DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                            .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                            .map(w -> getMarginalCDFDifference(instance.getWeights(), s.item, s.previousItemWeight, w)*immediateValueFunction.apply(s, x, w)+ 
                                                                      ((s.item < instance.getItems() - 1) ? getMarginalCDFDifference(instance.getWeights(), s.item, s.previousItemWeight, w)*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                            .sum())
                           .max()
                           .getAsDouble();
         double bestAction = Arrays.stream(s.getFeasibleActions())
                                .filter(x -> (
                                      (x == 0) ? immediateValueFunction.apply(s, 0.0, 0.0) + ((s.item < instance.getItems() - 1) ? f(stateTransition.apply(s, 0.0)) : 0) :
                                                 DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                             .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                             .map(w -> getMarginalCDFDifference(instance.getWeights(), s.item, s.previousItemWeight, w)*immediateValueFunction.apply(s, x, w)+ 
                                                                       ((s.item < instance.getItems() - 1) ? getMarginalCDFDifference(instance.getWeights(), s.item, s.previousItemWeight, w)*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                             .sum()) == val)
                                .findAny()
                                .getAsDouble();
         cacheActions.putIfAbsent(s, bestAction);
         return val;
      });
   }
   
   public static void main(String args[]) {

      SKPMultinormal instance = SKPMultinormal.getTestInstanceSpecialStructure();
      
      double truncationQuantile = 0.999;
      
      DSKPMultinormal dskp = new DSKPMultinormal(instance, truncationQuantile);
      
      System.out.println(GSONUtility.<DSKPMultinormalSolvedInstance>printInstanceAsJSON(dskp.solve()));
   }
   
   /*public static void main(String args[]) {
      SimpsonIntegrator simpson = new SimpsonIntegrator();
      double[] mu = {20, 30};
      double[][] sigma = {{16., 21.6}, {21.6, 36.}};
      MultiNormalDist dist = new MultiNormalDist(mu, sigma);

      UnivariateFunction uf = new MarginalCDF(dist, 20);

      double j = simpson.integrate(10000, uf, 30, 31);
      System.out.println("Trapezoid integral : " + j);
      System.out.println("Density : " + dist.density(new double[] {20,30}));
   }*/
   
   public static double getMarginalCDFDifference(MultiNormalDist dist, int coordinate, double realisedDemand, double value) {
      if(coordinate == 0) {
         NormalDist normal = new NormalDist(dist.getMu(0), Math.sqrt(dist.getSigma()[0][0]));
         return normal.cdf(value + d) - normal.cdf(value - d);
      } else {
         NormalDist normal = new NormalDist(dist.getMu(coordinate-1), Math.sqrt(dist.getSigma()[coordinate-1][coordinate-1]));
         UnivariateFunction pdf = new MarginalPDF(getMarginalDistribution(dist, new int[] {coordinate-1,coordinate}), realisedDemand);
         SimpsonIntegrator simpson = new SimpsonIntegrator();
         return simpson.integrate(100000, pdf, value - d, value + d)/(normal.cdf(realisedDemand + d) - normal.cdf(realisedDemand - d));
      }
   }
   
   public static MultiNormalDist getMarginalDistribution(MultiNormalDist dist, int[] coordinates) {
      double[] mu = new double[coordinates.length];
      for(int i = 0; i < coordinates.length; i++)
         mu[i] = dist.getMu(coordinates[i]);
      double[][] sigma = new double[coordinates.length][];
      for(int i = 0; i < coordinates.length; i++) {
         sigma[i] = new double[coordinates.length];
         for(int j = 0; j < coordinates.length; j++)
            sigma[i][j] = dist.getSigma()[coordinates[i]][coordinates[j]];
      }
      return new MultiNormalDist(mu, sigma);
   }
}

class MarginalPDF implements UnivariateFunction{
   MultiNormalDist dist;
   double realisedDemand;
   
   public MarginalPDF(MultiNormalDist dist, double realisedDemand) {
      this.dist = dist;
      this.realisedDemand = realisedDemand;
   }
   
   public double value(double x) {
      return dist.density(new double[]{realisedDemand, x});
   }
}


