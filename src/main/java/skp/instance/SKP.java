package skp.instance;

import java.util.Arrays;

import skp.utilities.hash.SHA;

public abstract class SKP {

   protected String instanceID;
   protected double[] expectedValues;
   protected double[] expectedWeights;
   protected double shortageCost;
   
   SKP(double[] expectedValues, double[] expectedWeights, double shortageCost){
      this.expectedValues = expectedValues;
      this.expectedWeights = expectedWeights;
      this.shortageCost = shortageCost;
   }
   
   void generateInstanceID() {
      String intHash = this.getClass().getSimpleName()+this.hashCode();
      instanceID = SHA.generateSHA256(intHash);
   }

   public String getInstanceID() {
      return this.instanceID;
   }

   public int getItems() {
      return this.expectedValues.length;
   }

   public double getShortageCost() {
      return this.shortageCost;
   }

   public double[] getExpectedValues() {
      return this.expectedValues;
   }
   
   /**
    * Returns a list of item indexes sorted by expected value to expected weight ratio in descending order.
    */
   public int[] getItemsOrderedByExpectedValueOverExpectedWeightRatio() {
      int n = expectedValues.length;
      Integer[] indexes = new Integer[n];
      for (int i = 0; i < n; i++) {
          indexes[i] = i;
      }

      Arrays.sort(indexes, (i, j) -> {
          double ratioI = expectedValues[i] / expectedWeights[i];
          double ratioJ = expectedValues[j] / expectedWeights[j];
          return Double.compare(ratioJ, ratioI); // descending order
      });

      int[] result = new int[n];
      for (int i = 0; i < n; i++) {
          result[i] = indexes[i];
      }

      return result;
  }
}