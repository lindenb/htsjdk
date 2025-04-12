package htsjdk.tribble.bgen;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
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
import java.util.zip.*;

/**
 * see also : https://github.com/limix/bgen
 *    https://github.com/limix/bgen/blob/main/bgen-file-format.pdf
 */

public class BGenCodec extends BinaryFeatureCodec<BGenFeature> {
    private BGenHeader header;
    private boolean debug = true;
    
    public BGenCodec() {
        }
    
    /** convert a long to int, throws a tribble exception if the number is greater than Integer.MAX_VALUE */
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
        free_area = null;
        
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
            number_of_individuals_M = longToUnsignedInt(binaryCodec.readUInt());
            }
        final BGenFeature ctx=new BGenFeature();
        

        ctx.variant_id = readStringUInt16(binaryCodec);
        System.err.println("variant_id is : "+ctx.variant_id);
        //
        ctx.rs_id = readStringUInt16(binaryCodec);
        System.err.println("rs_id is : "+ctx.rs_id);
        //
        ctx.chrom = readStringUInt16(binaryCodec);
        System.err.println("contig is : "+ctx.chrom);
        
        ctx.position_int32 = longToUnsignedInt(binaryCodec.readUInt());
        
      
        
        //ctx.genotypes=new BGenGenotype[num_samples];
        
        switch(header.getLayout()) {
            case 1 : readLayout1(ctx, binaryCodec); break;
            case 2 : readLayout2(ctx, binaryCodec); break;
            default: throw new TribbleException("unknown layout type (" + header.getLayout() + ")");
            }
        return ctx;
        }
    
    /** fully read 'len' bytes */
    private byte[] readNBytes(final BinaryCodec binaryCodec, final int len) throws IOException {
    	if(len<0) throw new IllegalArgumentException("negative bytes length:" + len);
        final byte[] bytes = new byte[len];
        binaryCodec.readBytes(bytes);
        return bytes;
    	}
    
    /**
     * read a string of length 'len'
     */
    private String readString(final BinaryCodec binaryCodec, final int len) throws IOException {
        return new String(readNBytes(binaryCodec, len));
    	}
    	
    private String[] readNAlleles(final BinaryCodec binaryCodec, int num_alleles) throws IOException {
      final String[] alleles = new String[num_alleles];
      for(int k=0;k< num_alleles;++k) {
            alleles[k] = readStringUInt64(binaryCodec);
            System.err.println("alelele "+(k+1)+"/"+num_alleles+" is : "+alleles[k]);
            }
      return alleles;
      }     
            
    /** read the length of a string as unsigned long and read the string */
    private String readStringUInt64(final BinaryCodec binaryCodec) throws IOException {
    	return readString(binaryCodec, longToUnsignedInt(binaryCodec.readUInt()));
    	}
    /** read the length of a string as unsigned short and read the string */
    private String readStringUInt16(final BinaryCodec binaryCodec) throws IOException {
    	return readString(binaryCodec, binaryCodec.readUShort());
    	}
    
    private void readLayout1(final BGenFeature ctx, final BinaryCodec binaryCodec) throws IOException {
    	// for layout 1 the number of alleles is always 2
	ctx.alleles = readNAlleles(binaryCodec, 2);
	throw new IOException("cannot read layout 1");
    	}
    

    
    private void readLayout2(final BGenFeature ctx, final BinaryCodec binaryCodec) throws IOException {
    	  final int  num_alleles  = binaryCodec.readUShort(); 
    	  System.err.println(" num_alleles is : "+num_alleles);
          ctx.alleles = readNAlleles(binaryCodec, num_alleles);
	 // The total length C of the rest of the data for this variant. Seeking forward this many bytes takes you to the next variant data block. 
    	 final int C_data_length = longToUnsignedInt(binaryCodec.readUInt());
         System.err.println(" data_length is : "+C_data_length);
         System.err.println(" getCompressedSNPBlocks is : "+ this.header.getCompressedSNPBlocks());
         final int uncompressed_data_length;
         if(header.getCompressedSNPBlocks()==0) {
              uncompressed_data_length = C_data_length;
              } else {
               uncompressed_data_length = (int)binaryCodec.readUInt();
              }
         System.err.println(" unciompressed data_length is : "+uncompressed_data_length);
         ctx.compressed_data = readNBytes(binaryCodec, C_data_length-4);
	
            try(InputStream zin = uncompress(new ByteArrayInputStream(ctx.compressed_data),header.getCompressedSNPBlocks())) {
            	int n_reads = 0;
                final BinaryCodec c2 = new BinaryCodec(zin);
                ctx.n_samples = longToUnsignedInt(c2.readUInt());
                n_reads += 4;
                if(ctx.n_samples <0) throw new TribbleException("ctx.n_samples <0 ("+ctx.n_samples +")");
                System.err.println(" :n_samples2 is : "+ctx.n_samples);
                
                
                final int n_alleles = (int)c2.readUShort();
                System.err.println(" :n_alleles is : "+n_alleles);
                ctx.min_ploidy  = (int)c2.readUByte();
                if(ctx.min_ploidy  <0 || ctx.min_ploidy >63) throw new TribbleException("bad ctx.min_ploidy  ("+ctx.min_ploidy +")");
                System.err.println(" :min_ploidy is : "+ ctx.min_ploidy);
                n_reads++;
                
                ctx.max_ploidy  = (int)c2.readUByte();
                System.err.println(" :max_ploidy is : "+ctx.max_ploidy);
                if(ctx.max_ploidy  <0 || ctx.max_ploidy >63) throw new TribbleException("bad ctx.max_ploidy  ("+ctx.max_ploidy +")");
                 n_reads++;
                
                byte[] ploidy_array = readNBytes(c2, ctx.n_samples);
                for(int x=0;x < ploidy_array.length ;x++) {System.err.print("ploidy["+(x+1)+"]="+(int)ploidy_array[x]+";");} System.err.println();
                n_reads += ctx.n_samples;
                
                boolean phased  = c2.readUByte()==1;
                System.err.println(" :phased is : "+phased);
                n_reads++;
                
                final int B=  (int)c2.readUByte();
                System.err.println(" :B is : "+B+" so number of bits is "+(B*ctx.n_samples));
                if(B  <0 || B > 32) throw new TribbleException("bad num of bits  ("+ B +")");
                n_reads++;
                
                
               
                
                System.err.println("bits read n="+(B*ctx.n_samples));
                
                //byte[] remains= readNBytes(c2, uncompressed_data_length - n_reads);
                int storage_gt = (int)Math.ceil((B*ctx.n_samples)/8.0);
                //byte[] storage= new byte[storage_gt];
                //c2.readBytes(storage);
                //BitSet genotypebits = BitSet.valueOf(storage);
                if(phased) {
                     int c;
		        BitReader br=new BitReader(c2.getInputStream());
		        for(int hap=0;hap<2;++hap) {
		        while((c=br.read())!=-1) {
		        	System.err.print(c);
		        	}
		        System.err.print("and then");
		        System.err.println();
		        }
                    }
                else
                    {
                    }
                c2.close();
           	}
    	}
    
    private static byte[] uncompress(byte[] array,int type) throws IOException {
    	try(InputStream in = uncompress(new ByteArrayInputStream(array),type)) {
        	try(ByteArrayOutputStream os = new ByteArrayOutputStream()) {
			htsjdk.samtools.util.IOUtil.copyStream(in, os);
			os.flush();
			return os.toByteArray();
			}
        	}
        }
        
    private static class BitReader {
    	private InputStream in;
    	byte curr;
    	int offset=8;
    	BitReader(InputStream in) {
    		this.in = in;
    		}
    	int read() throws IOException {
    		if(offset>=8) {
    			int c= in.read();
    			if(c==-1) return -1;
    			this.curr=(byte)c;
    			offset=0;
    			}
    		int bit = (curr >> (7 - offset)) & 1;
       	offset++;
        	return bit;
    		}
    		/*
    	    public int readNumber(int n) throws IOException {
		int result = 0;
		for (int i = 0; i < n; i++) {
			int bit = read();
			if (bit == -1) {
				return -1; // End of stream reached before n bits
			    }
			result = (result << 1) | bit;
			}
		  return result;
		  }*/
	public void close() throws IOException {
		in.close();
		}
    	}
    
    private static InputStream uncompress(InputStream is,int type) throws IOException {
        switch(type) {
            case 0: return is;
            case 1 : return new java.util.zip.InflaterInputStream(is);
            case 2 : return new ZstdCompressorInputStream(is);
            default: throw new IOException("unsupported compression type:"+type);
            }
        }
}
