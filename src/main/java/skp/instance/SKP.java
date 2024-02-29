package skp.instance;

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
}