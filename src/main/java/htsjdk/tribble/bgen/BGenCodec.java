package htsjdk.tribble.bgen;

import java.io.ByteArrayInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import org.apache.commons.compress.compressors.zstandard.ZstdCompressorInputStream;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.BinaryFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.Allele;

/**
 * 
 * 
 * 
 * 
 * see also : https://github.com/limix/bgen
 *
 */

public class BGenCodec extends BinaryFeatureCodec<BGenFeature> {
    
    static class BGenLayout
    {
        int    layout;
        int    nsamples;
        int    nalleles;
        int     phased;
        int     nbits;
        uint8_t*    ploidy_missingness;
        int    ncombs;
        int     min_ploidy;
        int     max_ploidy;
        byte[]       chunk;
        byte[] chunk_ptr;
        long    offset;
        };

    
    private BGenHeader header;
    public BGenCodec() {
        }
    
    static int longToUnsignedInt(long n) {
        if(n>Integer.MAX_VALUE || n<0) throw new TribbleException("Cannot convert "+n+"L to int");
        return (int)n;
    }
    
    @Override
    public boolean canDecode(final String path) {
        return path.endsWith(FileExtensions.BGEN);
    }
    @Override
    public Feature decodeLoc(PositionalBufferedStream source) throws IOException {
        return decode(source);
    }
    
    @Override
    public FeatureCodecHeader readHeader(PositionalBufferedStream source) throws IOException {
        final BGenHeader header = readActualHeader(source);
        return new FeatureCodecHeader(header,source.getPosition());
        }
    

    public BGenHeader readActualHeader(final PositionalBufferedStream seekableStream) {
        if(this.header!=null) throw new TribbleException("readRealHeader was already called");
        this.header = new BGenHeader();
        @SuppressWarnings("resource")
        final BinaryCodec binaryCodec = new BinaryCodec(seekableStream);
        header.variants_data_block_offset =  binaryCodec.readInt();
        
        int header_block_size =  longToUnsignedInt(binaryCodec.readUInt());
        header.num_variants_blocks =  binaryCodec.readUInt();
        final int num_samples = binaryCodec.readInt();
        
        final byte magic[]=new byte[4];
        binaryCodec.readBytes(magic);
        if(!(Arrays.equals(magic, "bgen".getBytes()) ||Arrays.equals(magic, new byte[] {0,0,0,0}))) {
            throw new TribbleException("Cannot decode magic "+magic[0]+"/"+magic[1]+"/"+magic[2]+"/"+magic[3]);
            }
       
        //Free data area. This could be used to store, for example, identifying information about the file
        byte[] free_area = new byte[header_block_size-20];
        binaryCodec.readBytes(free_area);
        // A set of flags, with bits numbered as for an unsigned integer. See below for flag definitions.
        byte flags[]=new byte[4];
        binaryCodec.readBytes(flags);
        header.headerFlags = BitSet.valueOf(flags);
      
    
        
        final boolean sample_identifiers_flags =  header.headerFlags.get(31);
        if(sample_identifiers_flags) {
            /* final long sample_block_size = */ binaryCodec.readUInt();        
            final int  number_of_samples = longToUnsignedInt(binaryCodec.readUInt());
            
            final List<String> samples= new ArrayList<>(number_of_samples);
            for(int i=0;i< number_of_samples;++i) {
                final int sn_len = binaryCodec.readUShort();
                final byte[] sample_name_as_bytes = new byte[sn_len];
                binaryCodec.readBytes(sample_name_as_bytes);
                samples.add(new String(sample_name_as_bytes));
                }
            header.samples= Collections.unmodifiableList(samples);
            }
        return header;
        }
    
    @Override
    public Class<BGenFeature> getFeatureType() {
        return BGenFeature.class;
    }
    
    // https://github.com/limix/bgen/blob/main/src/layout2.c
    
    @Override
    public BGenFeature decode(PositionalBufferedStream source) throws IOException {
        if(this.header==null) throw new TribbleException("readRealHeader was not called");
        @SuppressWarnings("resource")
        final BinaryCodec binaryCodec = new BinaryCodec(source);

        
        
        int number_of_individuals_M = -1;
        
        if(this.header.getLayout()==1) {
            try {
                number_of_individuals_M = longToUnsignedInt(binaryCodec.readUInt());
                }
            catch(EOFException err) {
                return null;
                }
            }
        BGenFeature ctx=new BGenFeature();
        

        int id_len = binaryCodec.readUShort();
        byte[] variant_id_as_bytes = new byte[id_len];
        binaryCodec.readBytes(variant_id_as_bytes);
        ctx.variant_id = new String(variant_id_as_bytes);
        System.err.println("variant_id is : "+ctx.variant_id);
        //
        id_len = binaryCodec.readUShort();
        variant_id_as_bytes = new byte[id_len];
        binaryCodec.readBytes(variant_id_as_bytes);
        ctx.rs_id = new String(variant_id_as_bytes);
        System.err.println("rs_id is : "+ctx.rs_id);
        //
        id_len = binaryCodec.readUShort();
        variant_id_as_bytes = new byte[id_len];
        binaryCodec.readBytes(variant_id_as_bytes);
        ctx.chrom = new String(variant_id_as_bytes);
        System.err.println("contig is : "+ctx.chrom);
        
        ctx.position_int32 = longToUnsignedInt(binaryCodec.readUInt());
        int  num_alleles;
        if(this.header.getLayout()!=1) {
            num_alleles = binaryCodec.readUShort();
            }
        else
            {
            num_alleles=2;
            }
        ctx.alleles = new Allele[num_alleles];

        for(int k=0;k< num_alleles;++k) {
            int len_allele_str = longToUnsignedInt(binaryCodec.readUInt());
            variant_id_as_bytes = new byte[len_allele_str];
            binaryCodec.readBytes(variant_id_as_bytes);
            ctx.alleles[k]=Allele.create(variant_id_as_bytes, k==0);
            System.err.println("alelele "+(k+1)+"/"+num_alleles+" is : "+ctx.alleles[k]);
            }
        
        ctx.genotypes=new BGenGenotype[num_samples];
        
        if(header.getLayout()==1) {
            throw new IllegalArgumentException("cannot read layout 1");
            }
        else
            {
            int data_length = (int)binaryCodec.readUInt();
            System.err.println(" data_length is : "+data_length);
            int total_length_D=data_length;
            if(header.getCompressedSNPBlocks()!=0) {
                total_length_D = (int)binaryCodec.readUInt();
                }
            byte[] compressed_data = new byte[data_length-4];
            binaryCodec.readBytes(compressed_data);
            CompressedGenotypeData gt=new CompressedGenotypeData(compressed_data);
            
            
            try(InputStream zin = uncompress(new ByteArrayInputStream(compressed_data),header.getCompressedSNPBlocks())) {
                BinaryCodec c2 = new BinaryCodec(zin);
                ctx.n_samples = (int)c2.readUInt();
                System.err.println(" :n_samples2 is : "+ctx.n_samples);
                int n_alleles = (int)c2.readUShort();
                System.err.println(" :n_alleles is : "+n_alleles);
                int min_ploidy  = (int)c2.readUByte();
                System.err.println(" :min_ploidy is : "+min_ploidy);
                int max_ploidy  = (int)c2.readUByte();
                System.err.println(" :max_ploidy is : "+max_ploidy);
                byte[] ploidy_array=new byte[ctx.n_samples];
                c2.readBytes(ploidy_array);
                int phased  = (int)c2.readUByte();
                System.err.println(" :phased is : "+phased);
                int B=  (int)c2.readUByte();
                System.err.println(" :B is : "+B+" so number of bits is "+(B*ctx.n_samples));
                int storage_gt = (int)Math.ceil((B*ctx.n_samples)/8.0);
                byte[] storage= new byte[storage_gt];
                c2.readBytes(storage);
                BitSet genotypebits = BitSet.valueOf(storage);
                if(phased==1) {
                    }
                else
                    {
                    }
                c2.close();
                }
            catch(IOException err) {
                err.printStackTrace();
                return null;
                }
            }
        
        }
    
    private static InputStream uncompress(InputStream is,int type) throws IOException {
        switch(type) {
            case 0: return is;
            case 2 : return new ZstdCompressorInputStream(is);
            default: throw new IOException("unsupported compression type:"+type);
            }
        }
    
    static  int read_ploidy(byte ploidy_miss) { 
        return ploidy_miss & 127;
        }
    
    static void read_phased_genotype_32(BGenLayout genotype, float[] probs)     
    {                                                                                         
        int nbits = genotype.nbits;                                                     
        int nalleles = genotype.nalleles;                                               
        int  max_ploidy = genotype.max_ploidy;        
        int probs_idx=0;
        float   denom = (float)((((uint64_t)1 << nbits)) - 1);                              
                                                                                              
        int sample_start = 0;                                                            
        for (int j = 0; j < genotype.nsamples; ++j) {                                   
            int ploidy = read_ploidy(genotype.ploidy_missingness[j]);                   
            float*  pend = probs + max_ploidy * nalleles;                                    
                                                                                              
            if (read_missingness(genotype.ploidy_missingness[j]) != 0) {                     
                set_array_nan##32(probs, (size_t)(pend - probs));                           
                probs = pend;                                                                 
                sample_start += nbits * ploidy * (nalleles - 1);                              
                continue;                                                                     
            }                                                                                 
                                                                                              
            int haplo_start = 0;                                                         
            for (int i = 0; i < ploidy; ++i) {                                            
                                                                                              
                int uip_sum = 0;                                                         
                int allele_start = 0;                                                    
                for (int ii = 0; ii < nalleles - 1; ++ii) {                              
                                                                                              
                    int ui_prob = 0;                                                     
                    int offset = sample_start + haplo_start + allele_start;              
                                                                                              
                    for (int bi = 0; bi < nbits; ++bi) {                                  
                                                                                              
                        if (get_bit(genotype.chunk_ptr, bi + offset)) {                      
                            ui_prob |= ((int)1 << bi);                                   
                        }                                                                     
                    }                                                                         
                                                                                              
                    probs[probs_idx] = (float)ui_prob / denom;                                         
                    ++probs_idx;                                                                  
                    uip_sum += ui_prob;                                                       
                    allele_start += nbits;                                                    
                }                                                                             
                probs[probs_idx] = (denom - (float)uip_sum) / denom;                                   
                ++probs_idx;                                                                      
                haplo_start += nbits * (nalleles - 1);                                        
            }                                                                                 
            sample_start += nbits * ploidy * (nalleles - 1);       
            
            for(int i=probs_idx;i< probs.length;++i) {
                probs[i]=Float.NaN;
                }
            probs = pend;                                                                     
        }                                                                                     
    }


    
}