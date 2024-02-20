package skp.sim;

import umontreal.ssj.rng.MRG32k3aL;

public abstract class Simulate {
   static final long[] seed = {1,2,3,4,5,6};
   MRG32k3aL randGenerator = new MRG32k3aL();
   
   public Simulate() {
      this.randGenerator.setSeed(seed);
   }
}
