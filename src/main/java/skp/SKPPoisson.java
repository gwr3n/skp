package skp;

import umontreal.ssj.probdist.PoissonDist;

public class SKPPoisson {
   double[] expectedValuesPerUnit;
   PoissonDist[] weights;
   int capacity;
   double shortageCost;
   
   public SKPPoisson(double[] expectedValuesPerUnit, PoissonDist[] weights, int capacity, double shortageCost) {
      this.expectedValuesPerUnit = expectedValuesPerUnit;
      this.weights = weights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
   }
   
   public int getItems() {
      return this.expectedValuesPerUnit.length;
   }
   
   public int getCapacity() {
      return this.capacity;
   }
   
   public double getShortageCost() {
      return this.shortageCost;
   }
   
   public double[] getExpectedValuesPerUnit() {
      return this.expectedValuesPerUnit;
   }
   
   public PoissonDist[] getWeights() {
      return this.weights;
   }
}
