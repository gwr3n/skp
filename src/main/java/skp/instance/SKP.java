package skp.instance;

import skp.utililities.hash.SHA;

public abstract class SKP {

   protected String instanceID;
   protected double[] expectedValuesPerUnit;
   protected double[] expectedWeights;
   protected double shortageCost;
   
   SKP(double[] expectedValuesPerUnit, double[] expectedWeights, double shortageCost){
      this.expectedValuesPerUnit = expectedValuesPerUnit;
      this.expectedWeights = expectedWeights;
      this.shortageCost = shortageCost;
   }
   
   void generateInstanceID() {
      String intHash = ""+this.hashCode();
      instanceID = SHA.generateSHA256(intHash);
   }

   public String getInstanceID() {
      return this.instanceID;
   }

   public int getItems() {
      return this.expectedValuesPerUnit.length;
   }

   public double getShortageCost() {
      return this.shortageCost;
   }

   public double[] getExpectedValuesPerUnit() {
      return this.expectedValuesPerUnit;
   }
}