package skp;

import java.math.BigInteger;
import java.util.Arrays;

public class KP {
   String instanceID;
   double[] valuesPerUnit;
   double[] weights;
   double capacity;
   double shortageCost;
   
   public KP(double[] valuesPerUnit, double[] weights, double capacity, double shortageCost) {
      this.valuesPerUnit = valuesPerUnit;
      this.weights = weights;
      this.capacity = capacity;
      this.shortageCost = shortageCost;
      
      generateInstanceID();
   }
   
   private void generateInstanceID() {
      String intHash = ""+this.hashCode();
      BigInteger hashCode = new BigInteger(intHash);
      instanceID = hashCode.toString(16);
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
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(this.valuesPerUnit) +
            Arrays.hashCode(this.weights) +
            Double.hashCode(this.capacity) + 
            Double.hashCode(this.shortageCost);
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof KP) {
         KP o = (KP) obj;
         return Arrays.equals(this.valuesPerUnit, o.valuesPerUnit) &&
                Arrays.equals(this.weights, o.weights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
}
