package skp.batch;

import umontreal.ssj.rng.MRG32k3aL;

public abstract class SKPBatch {

   static final long[] seed = {12345,54321,21435,53412,54321,14235};
   protected static final MRG32k3aL randGenerator = new MRG32k3aL();
}