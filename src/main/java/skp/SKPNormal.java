package skp;

import umontreal.ssj.probdist.NormalDist;

public class SKPNormal {
   double[] expectedValues;
   NormalDist[] weights;
   double capacity;
   double shortageCost;
   
   public SKPNormal(double[] expectedValues, NormalDist[] weights, double capacity, double shortageCost) {
      this.expectedValues = expectedValues;
      this.weights = weights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
   }
   
   public int getItems() {
      return this.expectedValues.length;
   }
   
   public double getCapacity() {
      return this.capacity;
   }
   
   public double getShortageCost() {
      return this.shortageCost;
   }
   
   public double[] getExpectedValues() {
      return this.expectedValues;
   }
   
   public NormalDist[] getWeights() {
      return this.weights;
   }
}
