package skp.instance;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;

import umontreal.ssj.probdist.BinomialDist;
import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.GammaDist;
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
      Distribution[] weights = {new NormalDist(44,10),
                                new PoissonDist(42),
                                new BinomialDist(73,0.5),
                                new GammaDist(15, 20),
                                new PoissonDist(71),
                                new PoissonDist(12),
                                new PoissonDist(13),
                                new PoissonDist(14),
                                new PoissonDist(23),
                                new PoissonDist(15)};
      int capacity = 100;
      int shortageCost = 10;
      return new SKPGenericDistribution(expectedValues, weights, capacity, shortageCost);
   }
   
   public static SKPGenericDistribution getTestInstanceLarge() {
      int objects = 35;
      Random rnd = new Random(1122);
      double[] expectedValues = IntStream.iterate(0, i -> i + 1)
                                         .limit(objects)
                                         .mapToDouble(i -> 100*rnd.nextDouble())
                                         .toArray();
      Distribution[] weights = IntStream.iterate(0, i -> i + 1)
                                        .limit(objects)
                                        .mapToObj(i -> new NormalDist(100*rnd.nextDouble(),10*rnd.nextDouble()))
                                        .toArray(NormalDist[]::new);
      int capacity = 10*50;
      int shortageCost = 10;
      return new SKPGenericDistribution(expectedValues, weights, capacity, shortageCost);
   }
}
