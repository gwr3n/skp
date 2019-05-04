package skp;

public class KP {
   double[] valuesPerUnit;
   double[] weights;
   double capacity;
   double shortageCost;
   
   public KP(double[] valuesPerUnit, double[] weights, double capacity, double shortageCost) {
      this.valuesPerUnit = valuesPerUnit;
      this.weights = weights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
   }
   
   public int getItems() {
      return this.valuesPerUnit.length;
   }
   
   public double getCapacity() {
      return this.capacity;
   }
   
   public double getShortageCost() {
      return this.shortageCost;
   }
   
   public double[] getValues() {
      return this.valuesPerUnit;
   }
   
   public double[] getWeights() {
      return this.weights;
   }
}
