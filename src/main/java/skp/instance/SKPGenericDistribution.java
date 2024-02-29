package skp.instance;

import java.util.Arrays;

import umontreal.ssj.probdist.BinomialDist;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

public class SKPGenericDistribution extends SKP{
   
   Distribution[] weights;
   double capacity;
   
   public SKPGenericDistribution(double[] expectedValues, Distribution[] weights, double capacity, double shortageCost) {
      super(expectedValues, Arrays.stream(weights).mapToDouble(w -> w.getMean()).toArray(), shortageCost);
      this.weights = weights;
      this.capacity = capacity;
      
      generateInstanceID();
   }

   public double getCapacity() {
      return this.capacity;
   }
   
   public Distribution[] getWeights() {
      return this.weights;
   }
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(new int[] {Arrays.hashCode(this.expectedValues),
                                        Arrays.hashCode(Arrays.stream(this.weights).map(w -> w.toString()).toArray(String[]::new)),
                                        Double.hashCode(this.capacity),
                                        Double.hashCode(this.shortageCost)});
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPGenericDistribution) {
         SKPGenericDistribution o = (SKPGenericDistribution) obj;
         return Arrays.equals(this.expectedValues, o.expectedValues) &&
                Arrays.equals(Arrays.stream(this.weights).map(w -> w.toString()).toArray(String[]::new), Arrays.stream(o.weights).map(w -> w.toString()).toArray(String[]::new)) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
   
   public static SKPGenericDistribution getTestInstance() {
      double[] expectedValues = {111, 111, 21, 117, 123, 34, 3, 121, 112, 12};
      Distribution[] expectedWeights = {new NormalDist(44,10),
                                        new PoissonDist(42),
                                        new BinomialDist(73,0.5),
                                        new PoissonDist(15),
                                        new PoissonDist(71),
                                        new PoissonDist(12),
                                        new PoissonDist(13),
                                        new PoissonDist(14),
                                        new PoissonDist(23),
                                        new PoissonDist(15)};
      int capacity = 100;
      int shortageCost = 10;
      return new SKPGenericDistribution(expectedValues, expectedWeights, capacity, shortageCost);
   }
}
