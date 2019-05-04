package skp;

import umontreal.ssj.probdistmulti.MultiNormalDist;

public class SKPMultiNormal {
   double[] expectedValues;
   MultiNormalDist weights;
   double capacity;
   double shortageCost;
   
   public SKPMultiNormal(double[] expectedValues, MultiNormalDist weights, double capacity, double shortageCost) {
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
   
   public MultiNormalDist getWeights() {
      return this.weights;
   }
   
}
