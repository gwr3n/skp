package skp;

import umontreal.ssj.probdist.PoissonDist;

public class SKPPoisson {
   double[] expectedValues;
   PoissonDist[] weights;
   int capacity;
   double shortageCost;
   
   public SKPPoisson(double[] expectedValues, PoissonDist[] weights, int capacity, double shortageCost) {
      this.expectedValues = expectedValues;
      this.weights = weights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
   }
   
   public int getItems() {
      return this.expectedValues.length;
   }
   
   public int getCapacity() {
      return this.capacity;
   }
   
   public double getShortageCost() {
      return this.shortageCost;
   }
   
   public double[] getExpectedValues() {
      return this.expectedValues;
   }
   
   public PoissonDist[] getWeights() {
      return this.weights;
   }
}
