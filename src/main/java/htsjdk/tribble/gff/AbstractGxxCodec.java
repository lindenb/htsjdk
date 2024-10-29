package htsjdk.tribble.gff;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.function.Predicate;
import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;

/**
 * Abstract Codec for parsing Gff3 files or GTF file
 */

abstract class AbstractGxxCodec extends AbstractFeatureCodec<Gff3Feature, LineIterator> {

    protected static final int NUM_FIELDS = 9;

    protected static final int CHROMOSOME_NAME_INDEX = 0;
    protected static final int ANNOTATION_SOURCE_INDEX = 1;
    protected static final int FEATURE_TYPE_INDEX = 2;
    protected static final int START_LOCATION_INDEX = 3;
    protected static final int END_LOCATION_INDEX = 4;
    protected static final int SCORE_INDEX = 5;
    protected static final int GENOMIC_STRAND_INDEX = 6;
    protected static final int GENOMIC_PHASE_INDEX = 7;
    protected static final int EXTRA_FIELDS_INDEX = 8;

    protected final DecodeDepth decodeDepth;

    /** filter to removing keys from the EXTRA_FIELDS column */
    protected final Predicate<String> filterOutAttribute;

    /**
     * @param decodeDepth a value from DecodeDepth
     * @param filterOutAttribute  filter to remove keys from the EXTRA_FIELDS column
     */
    protected AbstractGxxCodec(final DecodeDepth decodeDepth, final Predicate<String> filterOutAttribute) {
        super(Gff3Feature.class);
        this.decodeDepth = decodeDepth;
        this.filterOutAttribute = filterOutAttribute;
    }

    public enum DecodeDepth {
        DEEP ,
        SHALLOW
    }
    
    @Override
    public final Gff3Feature decode(final LineIterator lineIterator) throws IOException {
        return decode(lineIterator, decodeDepth);
    }

    protected abstract Gff3Feature decode(final LineIterator lineIterator, final DecodeDepth depth) throws IOException;

    @Override
    public final Feature decodeLoc(LineIterator lineIterator) throws IOException {
        return decode(lineIterator, DecodeDepth.SHALLOW);
    }
    
    static String extractSingleAttribute(final List<String> values) {
        if (values == null || values.isEmpty()) {
            return null;
        }

        if (values.size() != 1) {
            throw new TribbleException("Attribute has multiple values when only one expected");
        }
        return values.get(0);
    }

    @Override
    public final LineIterator makeSourceFromStream(final InputStream bufferedInputStream) {
        return new LineIteratorImpl(new SynchronousLineReader(bufferedInputStream));
    }

    @Override
    public final LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        return new AsciiLineReaderIterator(AsciiLineReader.from(bufferedInputStream));
    }

    @Override
    public final TabixFormat getTabixFormat() {
        return TabixFormat.GFF;
    }
}
