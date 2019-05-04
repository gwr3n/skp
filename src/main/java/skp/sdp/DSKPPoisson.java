package skp.sdp;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import umontreal.ssj.probdist.PoissonDist;

public class DSKPPoisson {
   
   int nbItems;
   PoissonDist[] distributions;
   int[] supportLB;
   int[] supportUB;

   public DSKPPoisson(int nbItems,
               PoissonDist[] distributions,
               int[] supportLB,
               int[] supportUB) {
      this.nbItems = nbItems;
      this.distributions = distributions;
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
                                                            .map(w -> distributions[s.item].prob((int)w)*immediateValueFunction.apply(s, x, w)+ 
                                                                      ((s.item < this.nbItems - 1) ? distributions[s.item].prob((int)w)*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                            .sum())
                           .max()
                           .getAsDouble();
         double bestAction = Arrays.stream(s.getFeasibleActions())
                                .filter(x -> (
                                      (x == 0) ? immediateValueFunction.apply(s, 0.0, 0.0) + ((s.item < this.nbItems - 1) ? f(stateTransition.apply(s, 0.0)) : 0) :
                                                 DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                             .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                             .map(w -> distributions[s.item].prob((int)w)*immediateValueFunction.apply(s, x, w)+ 
                                                                       ((s.item < this.nbItems - 1) ? distributions[s.item].prob((int)w)*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                             .sum()) == val)
                                .findAny()
                                .getAsDouble();
         cacheActions.putIfAbsent(s, bestAction);
         return val;
      });
   }

   public static void main(String args[]) {
      
      int nbItems = 5;
      double[] expectedValues = {100,100,100,100,100};
      double[] expectedWeights = {10,10,40,30,20};
      
      // Random variables
      
      double truncationQuantile = 0.999999999999999;

      PoissonDist[] distributions = IntStream.iterate(0, i -> i + 1)
                                              .limit(expectedWeights.length)
                                              //.mapToObj(i -> new NormalDist(expectedWeights[i],expectedWeights[i]*cv))
                                              .mapToObj(i -> new PoissonDist(expectedWeights[i]))
                                              .toArray(PoissonDist[]::new);
      int[] supportLB = IntStream.iterate(0, i -> i + 1)
                                    .limit(expectedWeights.length)
                                    .map(i -> (int)Math.round(distributions[i].inverseF(1-truncationQuantile)))
                                    .toArray();
      int[] supportUB = IntStream.iterate(0, i -> i + 1)
                                    .limit(expectedWeights.length)
                                    .map(i -> (int)Math.round(distributions[i].inverseF(truncationQuantile)))
                                    .toArray();
      
      int C = 100;
      int p = 100;
      
      DSKPPoisson dskp = new DSKPPoisson(nbItems, distributions, supportLB, supportUB);
      
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
            double value = action*expectedValues[state.item];
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
