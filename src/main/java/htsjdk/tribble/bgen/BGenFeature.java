package htsjdk.tribble.bgen;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;

public interface BGenFeature extends Feature {
   public int getNGenotypes();
   public BGenGenotype getGenotype(int index);
}
