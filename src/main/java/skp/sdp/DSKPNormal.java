package skp.sdp;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import umontreal.ssj.probdist.NormalDist;

public class DSKPNormal {
   
   int nbItems;
   NormalDist[] weightDistributions;
   int[] supportLB;
   int[] supportUB;
   double d = 0.5;
   
   public DSKPNormal(int nbItems,
                      NormalDist[] weightDistributions,
                      int[] supportLB,
                      int[] supportUB) {
      this.nbItems = nbItems;
      this.weightDistributions = weightDistributions;
      this.supportLB = supportLB;
      this.supportUB = supportUB;
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
                                                            .map(w -> (weightDistributions[s.item].cdf(w+d)-weightDistributions[s.item].cdf(w-d))*immediateValueFunction.apply(s, x, w)+ 
                                                                      ((s.item < this.nbItems - 1) ? (weightDistributions[s.item].cdf(w+d)-weightDistributions[s.item].cdf(w-d))*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                            .sum())
                           .max()
                           .getAsDouble();
         double bestAction = Arrays.stream(s.getFeasibleActions())
                                .filter(x -> (
                                      (x == 0) ? immediateValueFunction.apply(s, 0.0, 0.0) + ((s.item < this.nbItems - 1) ? f(stateTransition.apply(s, 0.0)) : 0) :
                                                 DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                             .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                             .map(w -> (weightDistributions[s.item].cdf(w+d)-weightDistributions[s.item].cdf(w-d))*immediateValueFunction.apply(s, x, w)+ 
                                                                       ((s.item < this.nbItems - 1) ? (weightDistributions[s.item].cdf(w+d)-weightDistributions[s.item].cdf(w-d))*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                             .sum()) == val)
                                .findAny()
                                .getAsDouble();
         cacheActions.putIfAbsent(s, bestAction);
         return val;
      });
   }
   
public static void main(String args[]) {
      
      int nbItems = 10;
      double[] expectedValuesPerUnit = {2.522727273, 2.642857143, 0.287671233, 7.8, 1.732394366, 2.833333333, 0.230769231, 8.642857143, 4.869565217, 0.8};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      
      // Random variables
      
      double truncationQuantile = 0.999999999999999;

      NormalDist[] weightDistributions = IntStream.iterate(0, i -> i + 1)
                                                  .limit(expectedWeights.length)
                                                  .mapToObj(i -> new NormalDist(expectedWeights[i],expectedWeights[i]*cv))
                                                  .toArray(NormalDist[]::new);
      int[] supportLB = IntStream.iterate(0, i -> i + 1)
                                    .limit(expectedWeights.length)
                                    .map(i -> (int)Math.round(weightDistributions[i].inverseF(1-truncationQuantile)))
                                    .toArray();
      int[] supportUB = IntStream.iterate(0, i -> i + 1)
                                    .limit(expectedWeights.length)
                                    .map(i -> (int)Math.round(weightDistributions[i].inverseF(truncationQuantile)))
                                    .toArray();
      
      int C = 100;
      int p = 100;
      
      DSKPNormal dskp = new DSKPNormal(nbItems, weightDistributions, supportLB, supportUB);
      
      /**
       * This function returns the set of actions associated with a given state
       */
      dskp.actionGenerator = state ->{
         return DoubleStream.iterate(0, num -> num + 1)
                            .limit(2)
                            .toArray();
      };
      
      /**
       * State transition function; given a state, an action and a random outcome, the function
       * returns the future state
       */
      dskp.stateTransition = (state, realisedWeight) -> 
         dskp.new State(state.item + 1, state.remainingCapacity - realisedWeight);
      
      /**
       * Immediate value function for a given state
       */
      dskp.immediateValueFunction = (state, action, realisedWeight) -> {
            double value = action*expectedValuesPerUnit[state.item]*realisedWeight;
            double cost = (state.item == nbItems - 1 ? p : 0)*Math.max(realisedWeight - state.remainingCapacity, 0);
            return value - cost;
         };
      
      /**
       * Initial problem conditions
       */
      int initialItem = 0;
      State initialState = dskp.new State(initialItem, C);
      
      /**
       * Run forward recursion and determine the probability of achieving the target wealth when
       * one follows an optimal policy
       */
      System.out.println("f_1("+C+")="+dskp.f(initialState));
      /**
       * Recover optimal action for item 1 when initial capacity at the beginning of period 1 is C.
       */
      System.out.println("b_1("+initialItem+")="+dskp.cacheActions.get(dskp.new State(initialItem, C)));
   }
}
