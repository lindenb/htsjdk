package htsjdk.tribble.gff;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Queue;
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

    protected final Queue<Gff3FeatureImpl> featuresToFlush = new ArrayDeque<>();
    protected final Queue<Gff3FeatureImpl> activeFeatures = new ArrayDeque<>();
    protected final Map<Integer, String> commentsWithLineNumbers = new LinkedHashMap<>();
    protected final DecodeDepth decodeDepth;

    /** filter to removing keys from the EXTRA_FIELDS column */
    protected final Predicate<String> filterOutAttribute;
    
    protected int currentLine = 0;

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

    /**
     * Gets map from line number to comment found on that line.  The text of the comment EXCLUDES the leading # which indicates a comment line.
     * @return Map from line number to comment found on line
     */
    public final Map<Integer, String> getCommentsWithLineNumbers() {
        return Collections.unmodifiableMap(new LinkedHashMap<>(commentsWithLineNumbers));
    }

    /**
     * Gets list of comments parsed by the codec.  Excludes leading # which indicates a comment line.
     * @return
     */
    public final List<String> getCommentTexts() {
        return Collections.unmodifiableList(new ArrayList<>(commentsWithLineNumbers.values()));
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
