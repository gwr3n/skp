package skp.instance;

import java.util.Arrays;

import umontreal.ssj.probdist.Distribution;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;

public class SKPGenericDistribution extends SKP{
   
   Distribution[] weights;
   double capacity;
   
   public SKPGenericDistribution(double[] expectedValuesPerUnit, Distribution[] weights, double capacity, double shortageCost) {
      super(expectedValuesPerUnit, Arrays.stream(weights).mapToDouble(w -> w.getMean()).toArray(), shortageCost);
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
      return Arrays.hashCode(new int[] {Arrays.hashCode(this.expectedValuesPerUnit),
                                        Arrays.hashCode(Arrays.stream(this.weights).map(w -> w.toString()).toArray(String[]::new)),
                                        Double.hashCode(this.capacity),
                                        Double.hashCode(this.shortageCost)});
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof SKPGenericDistribution) {
         SKPGenericDistribution o = (SKPGenericDistribution) obj;
         return Arrays.equals(this.expectedValuesPerUnit, o.expectedValuesPerUnit) &&
                Arrays.equals(Arrays.stream(this.weights).map(w -> w.toString()).toArray(String[]::new), Arrays.stream(o.weights).map(w -> w.toString()).toArray(String[]::new)) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
   
   public static SKPGenericDistribution getTestInstance() {
      //double[] expectedValuesPerUnit = {2.522727273, 2.642857143, 0.287671233, 7.8, 1.732394366, 2.833333333, 0.230769231, 8.642857143, 4.869565217, 0.8};
      //Distribution[] expectedWeights = {new PoissonDist(44),new PoissonDist(42),new PoissonDist(73),new PoissonDist(15),new PoissonDist(71),new PoissonDist(12),new PoissonDist(13),new PoissonDist(14),new PoissonDist(23),new PoissonDist(15)};
      double[] expectedValuesPerUnit = {2.522727273, 2.642857143, 0.287671233};
      Distribution[] expectedWeights = {new NormalDist(1,1),new NormalDist(1,1),new NormalDist(1,1)};
      int capacity = 50;
      int shortageCost = 100;
      return new SKPGenericDistribution(expectedValuesPerUnit, expectedWeights, capacity, shortageCost);
   }
}
