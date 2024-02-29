package skp.instance;

import java.util.Arrays;
import java.util.stream.IntStream;

import umontreal.ssj.probdist.PoissonDist;

public class SKPPoisson extends SKP {
   
   int capacity;
   
   public SKPPoisson(double[] expectedValues, double[] expectedWeights, int capacity, double shortageCost) {
      super(expectedValues, expectedWeights, shortageCost);
      this.capacity = capacity;
      
      generateInstanceID();
   }
   
   public SKPPoisson(double[] expectedValues, PoissonDist[] weights, int capacity, double shortageCost) {
      super(expectedValues, 
            IntStream.iterate(0, i -> i + 1).limit(weights.length)
                     .mapToDouble(i -> weights[i].getLambda())
                     .toArray(), 
            shortageCost);
      this.capacity = capacity;
      
      generateInstanceID();
   }
   
   public int getCapacity() {
      return this.capacity;
   }
   
   public PoissonDist[] getWeights() {
      return IntStream.iterate(0, i -> i + 1).limit(expectedWeights.length)
                      .mapToObj(i -> new PoissonDist(expectedWeights[i]))
                      .toArray(PoissonDist[]::new);
   }
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(new int[] {Arrays.hashCode(this.expectedValues),
                                        Arrays.hashCode(this.expectedWeights),
                                        Double.hashCode(this.capacity),
                                        Double.hashCode(this.shortageCost)});
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPPoisson) {
         SKPPoisson o = (SKPPoisson) obj;
         return Arrays.equals(this.expectedValues, o.expectedValues) &&
                Arrays.equals(this.expectedWeights, o.expectedWeights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
   
   public static SKPPoisson getTestInstance() {
      double[] expectedValues = {111, 111, 21, 117, 123, 34, 3, 121, 112, 12};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      int capacity = 100;
      int shortageCost = 10;
      return new SKPPoisson(expectedValues, expectedWeights, capacity, shortageCost);
   }
}
