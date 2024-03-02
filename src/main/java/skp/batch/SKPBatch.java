package skp.batch;

import umontreal.ssj.rng.MRG32k3aL;

public abstract class SKPBatch {

   static final long[] seed = {123456,654321,214365,653412,654321,142356};
   protected static final MRG32k3aL randGenerator = new MRG32k3aL();
}