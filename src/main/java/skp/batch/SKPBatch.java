package skp.batch;

import umontreal.ssj.rng.MRG32k3aL;

public abstract class SKPBatch {

   static final long[] seed = {123456, 246513, 264531, 425316, 635124, 326451};
   protected static final MRG32k3aL randGenerator = new MRG32k3aL();
}