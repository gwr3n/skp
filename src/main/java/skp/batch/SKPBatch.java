package skp.batch;

import umontreal.ssj.rng.MRG32k3aL;

public class SKPBatch {

   protected static final long[] seed = {1,2,3,4,5,6};
   protected static final MRG32k3aL randGenerator = new MRG32k3aL();

   public SKPBatch() {
      super();
   }

}