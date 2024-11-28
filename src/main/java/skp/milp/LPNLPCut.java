package skp.milp;

public class LPNLPCut {
   double[] X;  // Knapsack
   double LX;   // Value of loss function at X
   double[] dd; // Directional derivative
   
   public LPNLPCut(double[] X, double LX, double[] dd) {
      this.X = X;
      this.LX = LX;
      this.dd = dd;
   }
   
   public double[] getKnapsack() {
      return X;
   }
   
   public double getKnapsackLoss() {
      return this.LX;
   }
   
   public double[] getDirectionalDerivative() {
      return dd;
   }
}
