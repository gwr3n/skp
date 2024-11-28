package skp.milp;

public class Cut {
   int[] X;
   double rhs;
   
   public Cut(int[] X, double rhs) {
      this.X = X;
      this.rhs = rhs;
   }
   
   public int[] getKnapsack() {
      return X;
   }
   
   public double getRHS() {
      return rhs;
   }
}
