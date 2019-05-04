package skp;

import umontreal.ssj.probdist.NormalDist;

public class SKPNormal {
   double[] expectedValuesPerUnit;
   NormalDist[] weights;
   double capacity;
   double shortageCost;
   
   public SKPNormal(double[] expectedValuesPerUnit, NormalDist[] weights, double capacity, double shortageCost) {
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
   
   public NormalDist[] getWeights() {
      return this.weights;
   }
}
