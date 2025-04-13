package htsjdk.tribble.bgen;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.TribbleException;

public class BGenHeader {
    int variants_data_block_offset;
    List<String> samples;
    BitSet headerFlags;
    long num_variants_blocks;
    
    int getCompressedSNPBlocks() {
        return (this.headerFlags.get(1)?2:0)+
               (this.headerFlags.get(0)?1:0);
        }
    
    int getLayout() {
        return (this.headerFlags.get(5)?8:0)+
                (this.headerFlags.get(4)?4:0)+
                (this.headerFlags.get(3)?2:0)+
                (this.headerFlags.get(2)?1:0);
        }
    
    public List<String> getSamples() {
    	return this.samples;
    	}
    public int getNSamples() {
    	return this.samples.size();
    	}
   }
