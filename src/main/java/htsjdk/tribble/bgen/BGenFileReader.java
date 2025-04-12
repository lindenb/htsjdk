/**
 *
 *
 */
package htsjdk.tribble.bgen;



import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;

import java.io.File;
import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Iterator;


public class BGenFileReader  implements Closeable, Iterable<BGenFeature> {
    private final FeatureReader<BGenFeature> reader;
    public BGenFileReader(final Path path) {
    	this.reader = AbstractFeatureReader.getFeatureReader(
                path.toUri().toString(),
                new BGenCodec(),
                false);
    }
    
    
    /**
     * Returns an iterator over all records in this VCF/BCF file.
     */
    @Override
    public CloseableIterator<BGenFeature> iterator() {
        try {
            return reader.iterator();
        } catch (final IOException ioe) {
            throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
        }
    }

    
    @Override
    public void close() {
        try {
            this.reader.close();
        } catch (final IOException ioe) {
            throw new TribbleException("Could not close a variant context feature reader.", ioe);
        }
    }

}

