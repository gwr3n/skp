package skp.sim;

import skp.utilities.probability.SAMPLING;
import umontreal.ssj.rng.MRG32k3aL;

public abstract class Simulate {
   static final long[] seed = {12345,54321,21435,53412,54321,14235};
   MRG32k3aL randGenerator = new MRG32k3aL();
   
   public Simulate() {
      this.randGenerator.setSeed(seed);
   }
   
   public static final SAMPLING samplingStrategy = SAMPLING.LHS;
}
