package skp.sdp;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.sun.management.OperatingSystemMXBean;

import skp.SKPPoisson;
import skp.gson.GSONUtility;
import umontreal.ssj.probdist.PoissonDist;

@SuppressWarnings("restriction")
public class DSKPPoisson {
   
   int nbItems;
   SKPPoisson instance;
   PoissonDist[] weightDistributions;
   int[] supportLB;
   int[] supportUB;
   
   protected OperatingSystemMXBean osMBean;
   protected long nanoBefore;
   protected long nanoAfter;

   public DSKPPoisson(SKPPoisson instance, double truncationQuantile) {
      this.instance = instance;
      this.nbItems = instance.getItems();
      this.weightDistributions = instance.getWeights();
      supportLB = IntStream.iterate(0, i -> i + 1)
                           .limit(instance.getWeights().length)
                           .map(i -> (int)Math.round(weightDistributions[i].inverseF(1-truncationQuantile)))
                           .toArray();
      supportUB = IntStream.iterate(0, i -> i + 1)
                           .limit(instance.getWeights().length)
                           .map(i -> (int)Math.round(weightDistributions[i].inverseF(truncationQuantile)))
                           .toArray();
      
      initialiseFunctionalInterfaces();
   }
   
   public void initialiseFunctionalInterfaces() {
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
      this.new State(state.item + 1, state.remainingCapacity - realisedWeight);
      
      /**
       * Immediate value function for a given state
       */
      this.immediateValueFunction = (state, action, realisedWeight) -> {
            double value = action*instance.getExpectedValuesPerUnit()[state.item]*realisedWeight;
            double cost = (state.item == nbItems - 1 ? instance.getShortageCost() : 0)*Math.max(realisedWeight - state.remainingCapacity, 0);
            return value - cost;
         };
   }
   
   public DSKPPoissonSolvedInstance solve() {
      /**
       * Initial problem conditions
       */
      int initialItem = 0;
      State initialState = this.new State(initialItem, instance.getCapacity());
      
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
      
      return new DSKPPoissonSolvedInstance(instance, solutionValue, solutionTimeMs, statesExplored);
   }
   
   class State{
      int item;                           
      double remainingCapacity;            
      
      public State(int item, 
                   double remainingCapacity){
         this.item = item;
         this.remainingCapacity = remainingCapacity;
      }
      
      public double[] getFeasibleActions(){
         return actionGenerator.apply(this);
      }
      
      @Override
      public int hashCode(){
         String hash = "";
         hash = (hash + this.item) + "_" + this.remainingCapacity;
         return hash.hashCode();
      }
      
      @Override
      public boolean equals(Object o){
         if(o instanceof State)
            return  ((State) o).item == this.item &&
                    ((State) o).remainingCapacity == this.remainingCapacity;
         else
            return false;
      }

      @Override
      public String toString(){
         return this.item + " " + this.remainingCapacity;
      }
   }
   
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
   
   Map<State, Double> cacheActions = new HashMap<>();
   Map<State, Double> cacheValueFunction = new HashMap<>();
   double f(State state){
      return cacheValueFunction.computeIfAbsent(state, s -> {
         double val= Arrays.stream(s.getFeasibleActions())
                           .map(x -> (x == 0) ? immediateValueFunction.apply(s, 0.0, 0.0) + ((s.item < this.nbItems - 1) ? f(stateTransition.apply(s, 0.0)) : 0) :
                                                DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                            .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                            .map(w -> weightDistributions[s.item].prob((int)w)*immediateValueFunction.apply(s, x, w)+ 
                                                                      ((s.item < this.nbItems - 1) ? weightDistributions[s.item].prob((int)w)*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                            .sum())
                           .max()
                           .getAsDouble();
         double bestAction = Arrays.stream(s.getFeasibleActions())
                                .filter(x -> (
                                      (x == 0) ? immediateValueFunction.apply(s, 0.0, 0.0) + ((s.item < this.nbItems - 1) ? f(stateTransition.apply(s, 0.0)) : 0) :
                                                 DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                             .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                             .map(w -> weightDistributions[s.item].prob((int)w)*immediateValueFunction.apply(s, x, w)+ 
                                                                       ((s.item < this.nbItems - 1) ? weightDistributions[s.item].prob((int)w)*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                             .sum()) == val)
                                .findAny()
                                .getAsDouble();
         cacheActions.putIfAbsent(s, bestAction);
         return val;
      });
   }

   public static void main(String args[]) {
      double[] expectedValuesPerUnit = {2.522727273, 2.642857143, 0.287671233, 7.8, 1.732394366, 2.833333333, 0.230769231, 8.642857143, 4.869565217, 0.8};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      int C = 100;
      int p = 100;

      SKPPoisson instance = new SKPPoisson(expectedValuesPerUnit, expectedWeights, C, p);
      
      double truncationQuantile = 0.999999999999999;
      
      DSKPPoisson dskp = new DSKPPoisson(instance, truncationQuantile);
      
      System.out.println(GSONUtility.<DSKPPoissonSolvedInstance>printInstanceAsGSON(dskp.solve()));
   }
}
