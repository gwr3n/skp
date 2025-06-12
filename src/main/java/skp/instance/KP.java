package skp.instance;

import java.util.Arrays;

import skp.utilities.hash.SHA;

public class KP {
   String instanceID;
   double[] values;
   double[] weights;
   double capacity;
   double shortageCost;
   
   public KP(double[] values, double[] weights, double capacity, double shortageCost) {
      this.values = values;
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
   
   @Override
   public int hashCode() {
      return Arrays.hashCode(new int[] {Arrays.hashCode(this.values),
                                        Arrays.hashCode(this.weights),
                                        Double.hashCode(this.capacity),
                                        Double.hashCode(this.shortageCost)});
   }
   
   @Override
   public boolean equals(Object obj) {
      if(obj instanceof KP) {
         KP o = (KP) obj;
         return Arrays.equals(this.values, o.values) &&
                Arrays.equals(this.weights, o.weights) &&
                this.capacity == o.capacity &&
                this.shortageCost == o.shortageCost;
      }else
         return false;
   }
   
   /**
    * Returns a list of item indexes sorted by expected value to expected weight ratio in descending order.
    */
   public int[] getItemsOrderedByValueOverWeightRatio() {
      int n = values.length;
      Integer[] indexes = new Integer[n];
      for (int i = 0; i < n; i++) {
          indexes[i] = i;
      }

      Arrays.sort(indexes, (i, j) -> {
          double ratioI = values[i] / weights[i];
          double ratioJ = values[j] / weights[j];
          return Double.compare(ratioJ, ratioI); // descending order
      });

      int[] result = new int[n];
      for (int i = 0; i < n; i++) {
          result[i] = indexes[i];
      }

      return result;
  }
}
