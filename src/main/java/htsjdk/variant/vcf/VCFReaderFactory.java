package htsjdk.variant.vcf;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.function.Supplier;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;

/**
 * A Factory creating instances of {@link VCFReader}.
 * @author Pierre Lindenbaum
 * 
 */
public class VCFReaderFactory {
    /** provides a user-defined VCFReaderFactory. This field is used in {@link #newInstance()}. */
    private static Supplier<VCFReaderFactory> DEFAULT_FACTORY_PROVIDER = null;

    /** VCF source should be indexed ? */
    private boolean requireIndex = true;
    
    /** set the user-defined VCFReaderFactory that is going to be used in {@link #newInstance()}.
     * Setting to <code>null</code> restore the default behavior. */
    public static void setDefaultProvider(final Supplier<VCFReaderFactory> provider) {
        DEFAULT_FACTORY_PROVIDER = provider;
    }
    
    /** protected contstructor */
    protected VCFReaderFactory() {
    }
    
    /** if 'true' the VCFReader will assert the presence of an index as desired */
    public VCFReaderFactory setRequireIndex(boolean requireIndex) {
        this.requireIndex = requireIndex;
        return this;
    }
    
    /** return whether the VCFReader should assert the presence of an index.
     * The default value should be 'true' */
    public boolean isRequireIndex() {
        return this.requireIndex;
    }
    
    /** creates a new instance of  {@link #VCFReaderFactory()}, using the following methods:
     * (1) The a non-null Provider is used if was set with {@link #setDefaultProvider(Supplier)};
     * (2) A new instance is created with the class defined in {@link #DEFAULT_FACTORY_PROVIDER} if it's not blank;
     * (3) The default implementation is used. This default implementation creates a new {@link VCFFileReader}.
     *  */
    public static VCFReaderFactory newInstance() {
        if (DEFAULT_FACTORY_PROVIDER != null) return DEFAULT_FACTORY_PROVIDER.get();
        
        if (!StringUtil.isBlank(Defaults.VCF_READER_FACTORY)) {
            try {
                return (VCFReaderFactory)Class.forName(Defaults.VCF_READER_FACTORY).newInstance();
            } catch (final Throwable error) {
                throw new RuntimeException("Cannot create a new instance of VCFReaderFactory("
                        + Defaults.VCF_READER_FACTORY + ").", error);
                }
            }
        return new VCFReaderFactory();
        }
        
    /**
     * Constructs a VCFReader that will or will not assert the presence of an index as desired.
     */
    public VCFReader open(final String pathOrUrl) {     
        if (IOUtil.isUrl(pathOrUrl)) {
            throw new IllegalArgumentException("VCFReaderFactory["
                + this.getClass().getName()
                + "] cannot create a VCFReader from an URL ("
                + pathOrUrl + ").");
            }
        return open(Paths.get(pathOrUrl));
        }
    
    /**
     * Constructs a VCFReader that will or will not assert the presence of an index as desired.
     */
    public VCFReader open(final File vcfFile) {
        return open(vcfFile.toPath());
    }
        
    /**
     * Constructs a VCFReader that will or will not assert the presence of an index as desired.
     */
    public VCFReader open(final Path vcfFile) {
        return new VCFFileReader(vcfFile, isRequireIndex());
    }
}
