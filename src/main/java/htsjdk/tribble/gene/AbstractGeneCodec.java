package htsjdk.tribble.gene;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.LocationAware;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.SequenceRegion;
import htsjdk.tribble.gtf.GtfBaseData;
import htsjdk.tribble.gtf.GtfCodec;
import htsjdk.tribble.gtf.GtfConstants;
import htsjdk.tribble.gtf.GtfFeature;
import htsjdk.tribble.gtf.GtfFeatureImpl;
import htsjdk.tribble.gtf.GtfWriter;
import htsjdk.tribble.gtf.GtfCodec.DecodeDepth;
import htsjdk.tribble.gtf.GtfCodec.GtfDirective;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.tribble.util.ParsingUtils;

public abstract class AbstractGeneCodec<FEATURE_TYPE extends Feature> extends AbstractFeatureCodec<FEATURE_TYPE, LineIterator> { 

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


    protected static final String IS_CIRCULAR_ATTRIBUTE_KEY = "Is_circular";

    protected static final String ARTEMIS_FASTA_MARKER = ">";

    private final Queue<GtfFeatureImpl> activeFeatures = new ArrayDeque<>();
    private final Queue<GtfFeatureImpl> featuresToFlush = new ArrayDeque<>();
    private final Map<String, Set<GtfFeatureImpl>> activeFeaturesWithIDs = new HashMap<>();
    private final Map<String, Set<GtfFeatureImpl>> activeParentIDs = new HashMap<>();

    private final Map<String, SequenceRegion> sequenceRegionMap = new LinkedHashMap<>();
    private final Map<Integer, String> commentsWithLineNumbers = new LinkedHashMap<>();

    private final static Log logger = Log.getInstance(GtfCodec.class);

    private boolean reachedFasta = false;

    private DecodeDepth decodeDepth;

    private int currentLine = 0;

    /** filter to removing keys from the EXTRA_FIELDS column */
    private final Predicate<String> filterOutAttribute;
    
    /**
     * @param decodeDepth a value from DecodeDepth
     * @param filterOutAttribute  filter to remove keys from the EXTRA_FIELDS column
     */
    protected AbstractGeneCodec(final Class<FEATURE_TYPE> clazz,final DecodeDepth decodeDepth, final Predicate<String> filterOutAttribute) {
        super(clazz);
        this.decodeDepth = decodeDepth;
        this.filterOutAttribute = filterOutAttribute;
        
    }

    public enum DecodeDepth {
        DEEP ,
        SHALLOW
    }
    
    @Override
    public GtfFeature decode(final LineIterator lineIterator) throws IOException {
        return decode(lineIterator, decodeDepth);
    }

    private GtfFeature decode(final LineIterator lineIterator, final DecodeDepth depth) throws IOException {
        currentLine++;
        /*
        Basic strategy: Load features into deque, create maps from a features ID to it, and from a features parents' IDs to it.  For each feature, link to parents using these maps.
        When reaching flush directive, fasta, or end of file, prepare to flush features by moving all active features to deque of features to flush, and clearing
        list of active features and both maps.  Always poll featuresToFlush to return any completed top level features.
         */
        if (!lineIterator.hasNext()) {
            //no more lines, flush whatever is active
            prepareToFlushFeatures();
            return featuresToFlush.poll();
        }

        final String line = lineIterator.next();

        if (reachedFasta) {
            //previously reached fasta, flush whatever is active
            prepareToFlushFeatures();
            return featuresToFlush.poll();
        }

        if (line.startsWith(ARTEMIS_FASTA_MARKER)) {
            //backwards compatability with Artemis is built into gff3 spec
            processDirective(GtfDirective.FASTA_DIRECTIVE, null);
            return featuresToFlush.poll();
        }

        if (line.startsWith(GtfConstants.COMMENT_START) && !line.startsWith(GtfConstants.DIRECTIVE_START)) {
            commentsWithLineNumbers.put(currentLine, line.substring(GtfConstants.COMMENT_START.length()));
            return featuresToFlush.poll();
        }

        if (line.startsWith(GtfConstants.DIRECTIVE_START)) {
            parseDirective(line);
            return featuresToFlush.poll();
        }



        final GtfFeatureImpl thisFeature = new GtfFeatureImpl(parseLine(line, currentLine, this.filterOutAttribute));
        activeFeatures.add(thisFeature);
        if (depth == DecodeDepth.DEEP) {
            //link to parents/children/co-features
            final List<String> parentIDs = thisFeature.getAttribute(GtfConstants.PARENT_ATTRIBUTE_KEY);
            final String id = thisFeature.getID();

            for (final String parentID : parentIDs) {

                final Set<GtfFeatureImpl> theseParents = activeFeaturesWithIDs.get(parentID);
                if (theseParents != null) {
                    for (final GtfFeatureImpl parent : theseParents) {
                        thisFeature.addParent(parent);
                    }
                }
                if (activeParentIDs.containsKey(parentID)) {
                    activeParentIDs.get(parentID).add(thisFeature);
                } else {
                    activeParentIDs.put(parentID, new HashSet<>(Collections.singleton(thisFeature)));
                }
            }

            if (id != null) {
                if (activeFeaturesWithIDs.containsKey(id)) {
                    for (final GtfFeatureImpl coFeature : activeFeaturesWithIDs.get(id)) {
                        thisFeature.addCoFeature(coFeature);
                    }
                    activeFeaturesWithIDs.get(id).add(thisFeature);
                } else {
                    activeFeaturesWithIDs.put(id, new HashSet<>(Collections.singleton(thisFeature)));
                }
            }

            if (activeParentIDs.containsKey(thisFeature.getID())) {
                for (final GtfFeatureImpl child : activeParentIDs.get(thisFeature.getID())) {
                    child.addParent(thisFeature);
                }
            }
        }

        validateFeature(thisFeature);
        if (depth == DecodeDepth.SHALLOW) {
            //flush all features immediatly
            prepareToFlushFeatures();
        }
        return featuresToFlush.poll();
    }


    /**
     * Parse attributes field for gff3 feature
     * @param attributesString attributes field string from line in gff3 file
     * @return map of keys to values for attributes of this feature
     * @throws UnsupportedEncodingException
     */
    static private Map<String, List<String>> parseAttributes(final String attributesString) throws UnsupportedEncodingException {
        if (attributesString.equals(GtfConstants.UNDEFINED_FIELD_VALUE)) {
            return Collections.emptyMap();
        }
        final Map<String, List<String>> attributes = new LinkedHashMap<>();
        final List<String> splitLine = ParsingUtils.split(attributesString,GtfConstants.ATTRIBUTE_DELIMITER);
        for(String attribute : splitLine) {
            final List<String> key_value = ParsingUtils.split(attribute,GtfConstants.KEY_VALUE_SEPARATOR);
            if (key_value.size() != 2) {
                throw new TribbleException("Attribute string " + attributesString + " is invalid");
            }
            attributes.put(URLDecoder.decode(key_value.get(0).trim(), "UTF-8"), decodeAttributeValue(key_value.get(1).trim()));
        }
        return attributes;
    }

    private static GtfBaseData parseLine(final String line, final int currentLine, final Predicate<String> filterOutAttribute) {
        final List<String> splitLine = ParsingUtils.split(line, GtfConstants.FIELD_DELIMITER);

        if (splitLine.size() != NUM_FIELDS) {
            throw new TribbleException("Found an invalid number of columns in the given Gtf file at line + " + currentLine + " - Given: " + splitLine.size() + " Expected: " + NUM_FIELDS + " : " + line);
        }

        try {
            final String contig = URLDecoder.decode(splitLine.get(CHROMOSOME_NAME_INDEX), "UTF-8");
            final String source = URLDecoder.decode(splitLine.get(ANNOTATION_SOURCE_INDEX), "UTF-8");
            final String type = URLDecoder.decode(splitLine.get(FEATURE_TYPE_INDEX), "UTF-8");
            final int start = Integer.parseInt(splitLine.get(START_LOCATION_INDEX));
            final int end = Integer.parseInt(splitLine.get(END_LOCATION_INDEX));
            final double score = splitLine.get(SCORE_INDEX).equals(GtfConstants.UNDEFINED_FIELD_VALUE) ? -1 : Double.parseDouble(splitLine.get(SCORE_INDEX));
            final int phase = splitLine.get(GENOMIC_PHASE_INDEX).equals(GtfConstants.UNDEFINED_FIELD_VALUE) ? -1 : Integer.parseInt(splitLine.get(GENOMIC_PHASE_INDEX));
            final Strand strand = Strand.decode(splitLine.get(GENOMIC_STRAND_INDEX));
            final Map<String, List<String>> attributes = parseAttributes(splitLine.get(EXTRA_FIELDS_INDEX));
            /* remove attibutes matching 'filterOutAttribute' */
            attributes.keySet().removeIf(filterOutAttribute);
            return new GtfBaseData(contig, source, type, start, end, score, strand, phase, attributes);
        } catch (final NumberFormatException ex ) {
            throw new TribbleException("Cannot read integer value for start/end position from line " + currentLine + ".  Line is: " + line, ex);
        } catch (final IOException ex) {
            throw new TribbleException("Cannot decode feature info from line " + currentLine + ".  Line is: " + line, ex);
        }
    }

    /**
     * Get list of sequence regions parsed by the codec.
     * @return list of sequence regions
     */
    public List<SequenceRegion> getSequenceRegions() {
        return Collections.unmodifiableList(new ArrayList<>(sequenceRegionMap.values()));
    }

    /**
     * Gets map from line number to comment found on that line.  The text of the comment EXCLUDES the leading # which indicates a comment line.
     * @return Map from line number to comment found on line
     */
    public Map<Integer, String> getCommentsWithLineNumbers() {
        return Collections.unmodifiableMap(new LinkedHashMap<>(commentsWithLineNumbers));
    }

    /**
     * Gets list of comments parsed by the codec.  Excludes leading # which indicates a comment line.
     * @return
     */
    public List<String> getCommentTexts() {
        return Collections.unmodifiableList(new ArrayList<>(commentsWithLineNumbers.values()));
    }

    /**
     * If sequence region of feature's contig has been specified with sequence region directive, validates that
     * feature's coordinates are within the specified sequence region.  TribbleException is thrown if invalid.
     * @param feature
     */
    private void validateFeature(final GtfFeature feature) {
        if (sequenceRegionMap.containsKey(feature.getContig())) {
            final SequenceRegion region = sequenceRegionMap.get(feature.getContig());
            if (feature.getStart() == region.getStart() && feature.getEnd() == region.getEnd()) {
                //landmark feature
                final boolean isCircular = Boolean.parseBoolean(extractSingleAttribute(feature.getAttribute(IS_CIRCULAR_ATTRIBUTE_KEY)));
                region.setCircular(isCircular);
            }
            if (region.isCircular()? !region.overlaps(feature) : !region.contains(feature)) {
                throw new TribbleException("feature at " + feature.getContig() + ":" + feature.getStart() + "-" + feature.getEnd() +
                        " not contained in specified sequence region (" + region.getContig() + ":" + region.getStart() + "-" + region.getEnd());
            }
        }
    }

    @Override
    public Feature decodeLoc(LineIterator lineIterator) throws IOException {
        return decode(lineIterator, DecodeDepth.SHALLOW);
    }

    @Override
    public boolean canDecode(final String inputFilePath) {
        boolean canDecode;
        try {
            // Simple file and name checks to start with:
            Path p = IOUtil.getPath(inputFilePath);
            canDecode = FileExtensions.GFF3.stream().anyMatch(fe -> p.toString().endsWith(fe));

            if (canDecode) {

                // Crack open the file and look at the top of it:
                final InputStream inputStream = IOUtil.hasGzipFileExtension(p)? new GZIPInputStream(Files.newInputStream(p)) : Files.newInputStream(p);

                try ( BufferedReader br = new BufferedReader(new InputStreamReader(inputStream)) ) {

                    String line = br.readLine();

                    // First line must be GFF version directive
                    if (GtfDirective.toDirective(line) != GtfDirective.VERSION3_DIRECTIVE) {
                        return false;
                    }
                    while (line.startsWith(GtfConstants.COMMENT_START)) {
                        line = br.readLine();
                        if ( line == null ) {
                            return false;
                        }
                    }

                    // make sure line conforms to gtf spec
                    final List<String> fields = ParsingUtils.split(line, GtfConstants.FIELD_DELIMITER);

                    canDecode &= fields.size() == NUM_FIELDS;

                    if (canDecode) {
                        // check that start and end fields are integers
                        try {
                            final int start = Integer.parseInt(fields.get(3));
                            final int end = Integer.parseInt(fields.get(4));
                        } catch (NumberFormatException | NullPointerException nfe) {
                            return false;
                        }

                        // check for strand

                        final String strand = fields.get(GENOMIC_STRAND_INDEX);
                        canDecode &= strand.equals(Strand.POSITIVE.toString()) ||
                                strand.equals(Strand.NEGATIVE.toString()) ||
                                strand.equals(Strand.NONE.toString()) ||
                                strand.equals("?");
                    }
                }

            }
        }
        catch (FileNotFoundException ex) {
            logger.error(inputFilePath + " not found.");
            return false;
        }
        catch (final IOException ex) {
            return false;
        }

        return canDecode;
    }

    static List<String> decodeAttributeValue(final String attributeValue) {
        //split on VALUE_DELIMITER, then decode
        final List<String> splitValues = ParsingUtils.split(attributeValue, GtfConstants.VALUE_DELIMITER);

        final List<String> decodedValues = new ArrayList<>();
        for (final String encodedValue : splitValues) {
            try {
                decodedValues.add(URLDecoder.decode(encodedValue.trim(), "UTF-8"));
            } catch (final UnsupportedEncodingException ex) {
                throw new TribbleException("Error decoding attribute " + encodedValue, ex);
            }
        }

        return decodedValues;
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
    public FeatureCodecHeader readHeader(LineIterator lineIterator) {

        List<String> header = new ArrayList<>();
        while(lineIterator.hasNext()) {
            String line = lineIterator.peek();
            if (line.startsWith(GtfConstants.COMMENT_START)) {
                header.add(line);
                lineIterator.next();
            } else {
                break;
            }
        }

        return new FeatureCodecHeader(header, FeatureCodecHeader.NO_HEADER_END);
    }

    /**
     * Parse a directive line from a gff3 file
     * @param directiveLine
     * @throws IOException
     */
    private void parseDirective(final String directiveLine) throws IOException {
        final GtfDirective directive = GtfDirective.toDirective(directiveLine);
        if (directive != null) {
            processDirective(directive, directive.decode(directiveLine));

        } else {
            logger.warn("ignoring directive " + directiveLine);
        }
    }

    /**
     * Process a gff3 directive
     * @param directive the gff3 directive, indicated by a specific directive line
     * @param decodedResult the decoding of the directive line by the given directive
     */
    private void processDirective(final GtfDirective directive, final Object decodedResult) {
        switch (directive) {
            case VERSION3_DIRECTIVE:
                break;

            case SEQUENCE_REGION_DIRECTIVE:
                final SequenceRegion newRegion = (SequenceRegion) decodedResult;
                if (sequenceRegionMap.containsKey(newRegion.getContig())) {
                    throw new TribbleException("directive for sequence-region " + newRegion.getContig() + " included more than once.");
                }
                sequenceRegionMap.put(newRegion.getContig(), newRegion);
                break;

            case FLUSH_DIRECTIVE:
                prepareToFlushFeatures();
                break;

            case FASTA_DIRECTIVE:
                reachedFasta = true;
                break;

            default:
                throw new IllegalArgumentException( "Directive " + directive + " has been added to GtfDirective, but is not being handled by GtfCodec::processDirective.  This is a BUG.");

        }
    }

    /**
     * move active top level features to featuresToFlush.  clear active features.
     */
    private void prepareToFlushFeatures() {
        featuresToFlush.addAll(activeFeatures);
        activeFeaturesWithIDs.clear();
        activeFeatures.clear();
        activeParentIDs.clear();
    }

    @Override
    public LineIterator makeSourceFromStream(final InputStream bufferedInputStream) {
        return new LineIteratorImpl(new SynchronousLineReader(bufferedInputStream));
    }

    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        return new AsciiLineReaderIterator(AsciiLineReader.from(bufferedInputStream));
    }

    @Override
    public boolean isDone(final LineIterator lineIterator) {
        return !lineIterator.hasNext() && activeFeatures.isEmpty() && featuresToFlush.isEmpty();
    }

    @Override
    public void close(final LineIterator lineIterator) {
        //cleanup resources
        featuresToFlush.clear();
        activeFeaturesWithIDs.clear();
        activeFeatures.clear();
        activeParentIDs.clear();
        CloserUtil.close(lineIterator);
    }


    

}
