package skp.batch;

import umontreal.ssj.rng.MRG32k3aL;

public abstract class SKPBatch {

   static final long[] seed = {12345, 24513, 24531, 42531, 35124, 32451};
   protected static final MRG32k3aL randGenerator = new MRG32k3aL();
}