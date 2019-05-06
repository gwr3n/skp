package skp.sdp;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.sun.management.OperatingSystemMXBean;

import skp.instance.SKPNormal;
import skp.sdp.instance.DSKPNormalSolvedInstance;
import skp.utililities.gson.GSONUtility;

@SuppressWarnings("restriction")
public class DSKPNormal extends DSKP{
   
   SKPNormal instance;
   double d = 0.5;
   
   public DSKPNormal(SKPNormal instance, double truncationQuantile) {
      this.instance = instance;
      supportLB = IntStream.iterate(0, i -> i + 1)
                           .limit(instance.getWeights().length)
                           .map(i -> (int)Math.round(instance.getWeights()[i].inverseF(1-truncationQuantile)))
                           .toArray();
      supportUB = IntStream.iterate(0, i -> i + 1)
                           .limit(instance.getWeights().length)
                           .map(i -> (int)Math.round(instance.getWeights()[i].inverseF(truncationQuantile)))
                           .toArray();
      
      initialiseFunctionalInterfaces(instance);
   }
   
   public DSKPNormalSolvedInstance solve() {
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
      
      return new DSKPNormalSolvedInstance(instance, solutionValue, solutionTimeMs, statesExplored);
   }
   
   double f(State state){
      return cacheValueFunction.computeIfAbsent(state, s -> {
         double val= Arrays.stream(s.getFeasibleActions())
                           .map(x -> (x == 0) ? immediateValueFunction.apply(s, 0.0, 0.0) + ((s.item < instance.getItems() - 1) ? f(stateTransition.apply(s, 0.0)) : 0) :
                                                DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                            .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                            .map(w -> (instance.getWeights()[s.item].cdf(w+d)-instance.getWeights()[s.item].cdf(w-d))*immediateValueFunction.apply(s, x, w)+ 
                                                                      ((s.item < instance.getItems() - 1) ? (instance.getWeights()[s.item].cdf(w+d)-instance.getWeights()[s.item].cdf(w-d))*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                            .sum())
                           .max()
                           .getAsDouble();
         double bestAction = Arrays.stream(s.getFeasibleActions())
                                .filter(x -> (
                                      (x == 0) ? immediateValueFunction.apply(s, 0.0, 0.0) + ((s.item < instance.getItems() - 1) ? f(stateTransition.apply(s, 0.0)) : 0) :
                                                 DoubleStream.iterate(supportLB[s.item], k -> k + 1)
                                                             .limit(supportUB[s.item]-supportLB[s.item]+1)
                                                             .map(w -> (instance.getWeights()[s.item].cdf(w+d)-instance.getWeights()[s.item].cdf(w-d))*immediateValueFunction.apply(s, x, w)+ 
                                                                       ((s.item < instance.getItems() - 1) ? (instance.getWeights()[s.item].cdf(w+d)-instance.getWeights()[s.item].cdf(w-d))*f(stateTransition.apply(s, x == 0 ? 0 : w)) : 0))
                                                             .sum()) == val)
                                .findAny()
                                .getAsDouble();
         cacheActions.putIfAbsent(s, bestAction);
         return val;
      });
   }
   
   public static void main(String args[]) {

      SKPNormal instance = SKPNormal.getTestInstance();
      
      double truncationQuantile = 0.999999999999999;
      
      DSKPNormal dskp = new DSKPNormal(instance, truncationQuantile);
      
      System.out.println(GSONUtility.<DSKPNormalSolvedInstance>printInstanceAsGSON(dskp.solve()));
   }
}
