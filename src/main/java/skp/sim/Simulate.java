package skp.sim;

import skp.utilities.probability.SAMPLING;
import umontreal.ssj.rng.MRG32k3aL;

public abstract class Simulate {
   static final long[] seed = {123456, 246513, 264531, 425316, 635124, 326451};
   MRG32k3aL randGenerator = new MRG32k3aL();
   
   public Simulate() {
      this.randGenerator.setSeed(seed);
   }
   
   public static final SAMPLING samplingStrategy = SAMPLING.LHS;
}
