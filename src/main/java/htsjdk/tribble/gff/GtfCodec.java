package htsjdk.tribble.gff;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;

import htsjdk.samtools.util.Log;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.*;
import htsjdk.tribble.util.ParsingUtils;



import java.io.*;
import java.net.URLDecoder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * Codec for parsing GTF files, as defined in https://github.com/The-Sequence-Ontology/Specifications/blob/31f62ad469b31769b43af42e0903448db1826925/gff3.md
 * Note that while spec states that all feature types must be defined in sequence ontology, this implementation makes no check on feature types, and allows any string as feature type
 *
 * Each feature line in the Gff3 file will be emitted as a separate feature.  Features linked together through the "Parent" attribute will be linked through {@link Gff3Feature#getParents()}, {@link Gff3Feature#getChildren()},
 * {@link Gff3Feature#getAncestors()}, {@link Gff3Feature#getDescendents()}, amd {@link Gff3Feature#flatten()}.  This linking is not guaranteed to be comprehensive when the file is read for only features overlapping a particular
 * region, using a tribble index.  In this case, a particular feature will only be linked to the subgroup of features it is linked to in the input file which overlap the given region.
 */

public class GtfCodec extends AbstractGxxCodec {

    private final Map<String, Set<Gff3FeatureImpl>> activeFeaturesWithIDs = new HashMap<>();
    private final Map<String, Set<Gff3FeatureImpl>> activeParentIDs = new HashMap<>();

    private final static Log logger = Log.getInstance(GtfCodec.class);

    
    public GtfCodec() {
        this(DecodeDepth.DEEP);
    }

    public GtfCodec(final DecodeDepth decodeDepth) {
        this(decodeDepth, KEY -> false);
    }

    /**
     * @param decodeDepth a value from DecodeDepth
     * @param filterOutAttribute  filter to remove keys from the EXTRA_FIELDS column
     */
    public GtfCodec(final DecodeDepth decodeDepth, final Predicate<String> filterOutAttribute) {
       super(decodeDepth, filterOutAttribute);
        /* check required keys are always kept */
        for (final String key : new String[] {GtfConstants.GENE_ID, GtfConstants.TRANSCRIPT_ID}) {
            if (filterOutAttribute.test(key)) {
                throw new IllegalArgumentException("Predicate should always accept " + key);
                }
        }
    }

    protected Gff3Feature decode(final LineIterator lineIterator, final DecodeDepth depth) throws IOException {
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

        if (line.startsWith(GtfConstants.COMMENT_START)) {
            commentsWithLineNumbers.put(currentLine, line.substring(GtfConstants.COMMENT_START.length()));
            return featuresToFlush.poll();
        }

        final Gff3FeatureImpl thisFeature = new Gff3FeatureImpl(parseLine(line, currentLine, this.filterOutAttribute));
        activeFeatures.add(thisFeature);
        if (depth == DecodeDepth.DEEP) {
            //link to parents/children/co-features
            final List<String> parentIDs = thisFeature.getAttribute(GtfConstants.PARENT_ATTRIBUTE_KEY);
            final String id = thisFeature.getID();

            for (final String parentID : parentIDs) {

                final Set<Gff3FeatureImpl> theseParents = activeFeaturesWithIDs.get(parentID);
                if (theseParents != null) {
                    for (final Gff3FeatureImpl parent : theseParents) {
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
                    for (final Gff3FeatureImpl coFeature : activeFeaturesWithIDs.get(id)) {
                        thisFeature.addCoFeature(coFeature);
                    }
                    activeFeaturesWithIDs.get(id).add(thisFeature);
                } else {
                    activeFeaturesWithIDs.put(id, new HashSet<>(Collections.singleton(thisFeature)));
                }
            }

            if (activeParentIDs.containsKey(thisFeature.getID())) {
                for (final Gff3FeatureImpl child : activeParentIDs.get(thisFeature.getID())) {
                    child.addParent(thisFeature);
                }
            }
        }

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

    private static Gff3BaseData parseLine(final String line, final int currentLine, final Predicate<String> filterOutAttribute) {
        final List<String> splitLine = ParsingUtils.split(line, GtfConstants.FIELD_DELIMITER);

        if (splitLine.size() != NUM_FIELDS) {
            throw new TribbleException("Found an invalid number of columns in the given GTF file at line + " + currentLine + " - Given: " + splitLine.size() + " Expected: " + NUM_FIELDS + " : " + line);
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
            return new Gff3BaseData(contig, source, type, start, end, score, strand, phase, attributes);
        } catch (final NumberFormatException ex ) {
            throw new TribbleException("Cannot read integer value for start/end position from line " + currentLine + ".  Line is: " + line, ex);
        } catch (final IOException ex) {
            throw new TribbleException("Cannot decode feature info from line " + currentLine + ".  Line is: " + line, ex);
        }
    }

    @Override
    public boolean canDecode(final String inputFilePath) {
        boolean canDecode;
        try {
            // Simple file and name checks to start with:
            Path p = IOUtil.getPath(inputFilePath);
            canDecode = FileExtensions.GTF.stream().anyMatch(fe -> p.toString().endsWith(fe));

            if (canDecode) {

                // Crack open the file and look at the top of it:
                try( InputStream inputStream = IOUtil.hasGzipFileExtension(p)? new GZIPInputStream(Files.newInputStream(p)) : Files.newInputStream(p)) {
                    try ( BufferedReader br = new BufferedReader(new InputStreamReader(inputStream)) ) {
    
                        String line = br.readLine();
    
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
                                Integer.parseInt(fields.get(3));
                                Integer.parseInt(fields.get(4));
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
     * move active top level features to featuresToFlush.  clear active features.
     */
    private void prepareToFlushFeatures() {
        featuresToFlush.addAll(activeFeatures);
        activeFeaturesWithIDs.clear();
        activeFeatures.clear();
        activeParentIDs.clear();
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