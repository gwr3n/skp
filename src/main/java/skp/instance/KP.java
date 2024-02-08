package skp.instance;

import java.util.Arrays;

import skp.utililities.hash.SHA;

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
      String intHash = this.getClass().getSimpleName()+this.hashCode();
      instanceID = SHA.generateSHA256(intHash);
   }
   
   public String getInstanceID() {
      return this.instanceID;
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
      return Arrays.hashCode(new int[] {Arrays.hashCode(this.valuesPerUnit),
                                        Arrays.hashCode(this.weights),
                                        Double.hashCode(this.capacity),
                                        Double.hashCode(this.shortageCost)});
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
