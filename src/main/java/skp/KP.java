package skp;

public class KP {
   double[] values;
   double[] weights;
   double capacity;
   double shortageCost;
   
   public KP(double[] values, double[] weights, double capacity, double shortageCost) {
      this.values = values;
      this.weights = weights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
   }
   
   public int getItems() {
      return this.values.length;
   }
   
   public double getCapacity() {
      return this.capacity;
   }
   
   public double getShortageCost() {
      return this.shortageCost;
   }
   
   public double[] getValues() {
      return this.values;
   }
   
   public double[] getWeights() {
      return this.weights;
   }
}
