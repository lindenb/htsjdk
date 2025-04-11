package htsjdk.tribble.bgen;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;

public class BGenFeature implements Feature {
    String variant_id;
    String rs_id;
    String chrom;
    int position_int32;
    Allele[] alleles;
    byte[] ploidy;
    //
    int n_samples;
    
    @Override
    public String getContig() {
        return chrom;
        }
    @Override
    public int getStart() {
        return position_int32;
        }
    @Override
    public int getEnd() {
        return position_int32 + (alleles.length==0?0: alleles[0].length()-1);
        }
    
}
