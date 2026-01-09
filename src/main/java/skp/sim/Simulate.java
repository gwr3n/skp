package skp.sim;

import skp.utilities.probability.SAMPLING;
import umontreal.ssj.rng.MRG32k3aL;

public abstract class Simulate {
   static final long[] seed = {12345, 24513, 24531, 42531, 35124, 32451};
   MRG32k3aL randGenerator = new MRG32k3aL();
   
   public Simulate() {
      this.randGenerator.setSeed(seed);
   }
   
   public static final SAMPLING samplingStrategy = SAMPLING.LHS;
}
