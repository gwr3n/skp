package skp.sdp;

import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.DoubleStream;

import com.sun.management.OperatingSystemMXBean;

import skp.instance.SKP;

@SuppressWarnings("restriction")
public abstract class DSKP {
   
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
         this.new State(state.item + 1, state.remainingCapacity - realisedWeight);
      
      /**
       * Immediate value function for a given state
       */
      this.immediateValueFunction = (state, action, realisedWeight) -> {
            double value = action*instance.getExpectedValuesPerUnit()[state.item]*realisedWeight;
            double cost = (state.item == instance.getItems() - 1 ? instance.getShortageCost() : 0)*Math.max(realisedWeight - state.remainingCapacity, 0);
            return value - cost;
         };
   }
   
   Map<State, Double> cacheActions = new HashMap<>();
   Map<State, Double> cacheValueFunction = new HashMap<>();
   
   abstract double f(State state); 
   
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
}
