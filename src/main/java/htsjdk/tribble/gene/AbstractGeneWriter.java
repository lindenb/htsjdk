package htsjdk.tribble.gene;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.TribbleException;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.IOException;
import java.io.OutputStream;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;


/**
 * A class to write out gff3 files.  Features are added using {@link #addFeature(Gff3Feature)}, directives using {@link #addDirective(Gff3Codec.Gff3Directive)},
 * and comments using {@link #addComment(String)}.  Note that the version 3 directive is automatically added at creation, so should not be added separately.
 */
public abstract class AbstractGeneWriter<T extends GeneFeature<T>> implements Closeable {

    private final OutputStream out;
    private final AbstractGeneConstants geneConstants;

    
    protected AbstractGeneWriter(AbstractGeneConstants geneConstants,final Path path) throws IOException {
    	this.geneConstants = geneConstants;
    	validatePath(path);
    	final OutputStream outputStream = IOUtil.hasGzipFileExtension(path)? new BlockCompressedOutputStream(path.toFile()) : Files.newOutputStream(path);
        out = new BufferedOutputStream(outputStream);
        //start with version directive
        initialize();
    }

    public AbstractGeneWriter(AbstractGeneConstants geneConstants,final OutputStream stream) {
    	this.geneConstants = geneConstants;
        out = stream;
        initialize();
    }

    protected abstract void validatePath(final Path path) throws TribbleException;
    
    private void initialize() {
        try {
            writeWithNewLine(Gff3Codec.Gff3Directive.VERSION3_DIRECTIVE.encode(version));
        } catch (final IOException ex) {
            throw new TribbleException("Error writing version directive", ex);
        }
    }

    private void writeWithNewLine(final String txt) throws IOException {
        out.write(txt.getBytes());
        out.write(AbstractGeneConstants.END_OF_LINE_CHARACTER);
    }

    private void tryToWrite(final String string) {
        try {
            out.write(string.getBytes());
        } catch (final IOException ex) {
            throw new TribbleException("Error writing out string " + string, ex);
        }
    }

    protected void writeFirstEightFields(final T feature) throws IOException {
        writeJoinedByDelimiter(AbstractGeneConstants.FIELD_DELIMITER, this::tryToWrite, Arrays.asList(
                escapeString(feature.getContig()),
                escapeString(feature.getSource()),
                escapeString(feature.getType()),
                Integer.toString(feature.getStart()),
                Integer.toString(feature.getEnd()),
                feature.getScore() < 0 ? AbstractGeneConstants.UNDEFINED_FIELD_VALUE : Double.toString(feature.getScore()),
                feature.getStrand().toString(),
                feature.getPhase() < 0 ? AbstractGeneConstants.UNDEFINED_FIELD_VALUE : Integer.toString(feature.getPhase())
                )
        );
    }

    void writeAttributes(final Map<String, List<String>> attributes) throws IOException {
        if (attributes.isEmpty()) {
            out.write(AbstractGeneConstants.UNDEFINED_FIELD_VALUE.getBytes());
        }

        writeJoinedByDelimiter(AbstractGeneConstants.ATTRIBUTE_DELIMITER, e ->  writeKeyValuePair(e.getKey(), e.getValue()), attributes.entrySet());
    }

    void writeKeyValuePair(final String key, final List<String> values) {
        try {
            tryToWrite(key);
            out.write(AbstractGeneConstants.KEY_VALUE_SEPARATOR);
            writeJoinedByDelimiter(AbstractGeneConstants.VALUE_DELIMITER, v -> tryToWrite(escapeString(v)), values);
        } catch (final IOException ex) {
            throw new TribbleException("error writing out key value pair " + key + " " + values);
        }
    }

    private <T> void writeJoinedByDelimiter(final char delimiter, final Consumer<T> consumer, final Collection<T> fields) throws IOException {
        boolean isNotFirstField = false;
        for (final T field : fields) {
            if (isNotFirstField) {
                out.write(delimiter);
            } else {
                isNotFirstField = true;
            }

            consumer.accept(field);
        }
    }

    /***
     * add a feature
     * @param feature the feature to be added
     * @throws IOException
     */
    public void addFeature(final T feature) throws IOException {
        writeFirstEightFields(feature);
        out.write(AbstractGeneConstants.FIELD_DELIMITER);
        writeAttributes(feature.getAttributes());
        out.write(AbstractGeneConstants.END_OF_LINE_CHARACTER);
    }

    /***
     * escape a String.
     * Default behavior is to call {@link #encodeString(String)}
     * @param s the string to be escaped
     * @return the escaped string
     */
    protected String escapeString(final String s) {
        return encodeString(s);
    }

    static String encodeString(final String s) {
        try {
            //URLEncoder.encode is hardcoded to change all spaces to +, but we want spaces left unchanged so have to do this
            //+ is escaped to %2B, so no loss of information
            return URLEncoder.encode(s, "UTF-8").replace("+", " ");
        } catch (final UnsupportedEncodingException ex) {
            throw new TribbleException("Encoding failure", ex);
        }
    }

    /**
     * Add a directive with an object
     * @param directive the directive to be added
     * @param object the object to be encoded with the directive
     * @throws IOException
     */
    public void addDirective(final Gff3Codec.Gff3Directive directive, final Object object) throws IOException {
        if (directive == Gff3Codec.Gff3Directive.VERSION3_DIRECTIVE) {
            throw new TribbleException("VERSION3_DIRECTIVE is automatically added and should not be added manually.");
        }
        writeWithNewLine(directive.encode(object));
    }


    /**
     * Add comment line
     * @param comment the comment line (not including leading #)
     * @throws IOException
     */
    public void addComment(final String comment) throws IOException {
        out.write(AbstractGeneConstants.COMMENT_START.getBytes());
        writeWithNewLine(comment);
    }

    @Override
    public void close() throws IOException {
        out.close();
    }
}