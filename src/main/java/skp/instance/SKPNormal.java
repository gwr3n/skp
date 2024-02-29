package skp.instance;

import java.util.Arrays;
import java.util.stream.IntStream;

import umontreal.ssj.probdist.NormalDist;

public class SKPNormal extends SKP {

   double[] stdWeights;
   double capacity;
   
   public SKPNormal(double[] expectedValues, double[] expectedWeights, double[] stdWeights, double capacity, double shortageCost) {
      super(expectedValues, expectedWeights, shortageCost);
      this.stdWeights = stdWeights;
      this.capacity = capacity;
      
      generateInstanceID();
   }
   
   public SKPNormal(double[] expectedValues, double[] expectedWeights, double coefficientOfVariation, double capacity, double shortageCost) {
      this(expectedValues, 
            expectedWeights,
            IntStream.iterate(0, i -> i + 1)
                     .limit(expectedWeights.length)
                     .mapToDouble(i -> expectedWeights[i]*coefficientOfVariation)
                     .toArray(),
            capacity, 
            shortageCost);
   }
   
   public SKPNormal(double[] expectedValues, NormalDist[] weights, double capacity, double shortageCost) {
      this(expectedValues,
            IntStream.iterate(0, i -> i + 1)
                     .limit(weights.length)
                     .mapToDouble(i -> weights[i].getMean())
                     .toArray(),
            IntStream.iterate(0, i -> i + 1)
                     .limit(weights.length)
                     .mapToDouble(i -> weights[i].getSigma())
                     .toArray(),
            capacity,
            shortageCost);
   }
   
   public double getCapacity() {
      return this.capacity;
   }
   
   public NormalDist[] getWeights() {
      return IntStream.iterate(0, i -> i + 1).limit(expectedWeights.length)
                      .mapToObj(i -> new NormalDist(expectedWeights[i], stdWeights[i]))
                      .toArray(NormalDist[]::new);
   }
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(new int[] {Arrays.hashCode(this.expectedValues),
                                        Arrays.hashCode(this.expectedWeights),
                                        Arrays.hashCode(this.stdWeights),
                                        Double.hashCode(this.capacity),
                                        Double.hashCode(this.shortageCost)});
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPNormal) {
         SKPNormal o = (SKPNormal) obj;
         return Arrays.equals(this.expectedValues, o.expectedValues) &&
                Arrays.equals(this.expectedWeights, o.expectedWeights) &&
                Arrays.equals(this.stdWeights, o.stdWeights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
   
   public static SKPNormal getTestInstance() {
      double[] expectedValues = {111, 111, 21, 117, 123, 34, 3, 121, 112, 12};
      double[] expectedWeights = {44,42,73,15,71,12,13,14,23,15};
      double cv = 0.2;
      int capacity = 100;
      int shortageCost = 10;
      return new SKPNormal(expectedValues, expectedWeights, cv, capacity, shortageCost);
   }
}
