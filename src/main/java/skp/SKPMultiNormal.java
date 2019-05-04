package skp;

import umontreal.ssj.probdistmulti.MultiNormalDist;

public class SKPMultiNormal {
   double[] expectedValuesPerUnit;
   MultiNormalDist weights;
   double capacity;
   double shortageCost;
   
   public SKPMultiNormal(double[] expectedValuesPerUnit, MultiNormalDist weights, double capacity, double shortageCost) {
      this.expectedValuesPerUnit = expectedValuesPerUnit;
      this.weights = weights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
   }
   
   public int getItems() {
      return this.expectedValuesPerUnit.length;
   }
   
   public double getCapacity() {
      return this.capacity;
   }
   
   public double getShortageCost() {
      return this.shortageCost;
   }
   
   public double[] getExpectedValuesPerUnit() {
      return this.expectedValuesPerUnit;
   }
   
   public MultiNormalDist getWeights() {
      return this.weights;
   }
   
}
