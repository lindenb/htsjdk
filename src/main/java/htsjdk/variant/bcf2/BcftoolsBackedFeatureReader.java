package htsjdk.variant.bcf2;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class BcftoolsBackedFeatureReader implements FeatureReader<VariantContext> {
    private static String BCFTOOLS_EXE = null;
    private final VCFHeader header;
    private final boolean hasIndex;
    BcftoolsBackedFeatureReader(String bcfpath) {
        final Path path = Paths.get(bcfpath);
        IOUtil.assertFileIsReadable(path);
        final List<String> cmd = getCommandViewPreamble();
        cmd.add("--header-only");
        cmd.add(bcfpath);
       
        Process proc = null;
        InputStream in = null;
        try {
            proc =new ProcessBuilder(cmd).
                redirectError(Redirect.INHERIT).
                start();
            in  = proc.getInputStream();
            final VCFCodec codec = new VCFCodec();
            final LineIterator lineIter = codec.makeSourceFromStream(in);
            this.header = (VCFHeader)codec.readActualHeader(lineIter);
            if(proc.waitFor()!=0) throw new RuntimeIOException("Cannot read header from "+bcfpath);
            }   
        catch(final Throwable err) {
            throw new RuntimeIOException(err);
            }
        finally
            {
            CloserUtil.close(in);
            proc.destroy();
            }
        
        final Path tbiIdx =Paths.get(bcfpath+FileExtensions.TABIX_INDEX) ;
        final Path csiIdx =Paths.get(bcfpath+FileExtensions.CSI) ;
        
        this.hasIndex = (Files.exists(tbiIdx) || Files.exists(csiIdx));
        
    }
    
    @Override
    public boolean isQueryable() {
        return hasIndex;
    }

    private static String getBcfToolsExecutable() {
        if(!StringUtil.isBlank(BCFTOOLS_EXE)) return BCFTOOLS_EXE;
        
        BCFTOOLS_EXE =  System.getProperty("bcftools.exe");
        
        if(StringUtil.isBlank(BCFTOOLS_EXE)) {
            final String linux_path = System.getenv("PATH");
            if(!StringUtil.isBlank(linux_path)) {
                for(final String dir: linux_path.split(Pattern.quote(File.pathSeparator))) {
                    if(StringUtil.isBlank(dir)) continue;
                    final File exe = new File(dir,"bcftools");
                    if(!exe.exists()) continue;
                    if(exe.isDirectory()) continue;
                    if(!exe.canExecute()) continue;
                    BCFTOOLS_EXE = exe.getPath();
                    break;
                    }
                }
            }
        
        if(StringUtil.isBlank(BCFTOOLS_EXE)) {
            BCFTOOLS_EXE = "bcftools";
            }
        return BCFTOOLS_EXE;
        }
    
    @Override
    public Object getHeader() {
        return this.header;
    }
    
    private static List<String> getCommandViewPreamble() {
        final List<String> cmd = new ArrayList<>();
        cmd.add(BcftoolsBackedFeatureReader.getBcfToolsExecutable());
        cmd.add("view");
        cmd.add("-O");
        cmd.add("v");
        return cmd;
    }

    private VCFIterator commandToVcfIterator(final List<String> cmd,final VCFHeader header)  {
        Process proc = null;
        try {
            proc =new ProcessBuilder(cmd).
                redirectError(Redirect.INHERIT).
                start();
            }
        catch(final Throwable err) {
            throw new RuntimeIOException(err);
            }
        VCFIterator delegate = null;
        InputStream in = null;
        try {
            in  = proc.getInputStream();
            delegate = new  VCFIteratorBuilder().open(in);
            }
        catch(final Throwable err) {
            CloserUtil.close(in);
            proc.destroy();
            throw new RuntimeIOException(err);
            }
        final BcfToolsVcfIterator biter = new BcfToolsVcfIterator(VCFHeader header,delegate,proc);
        return biter;
        }

    
    
    private class BcfToolsViewImpl  {
        private final String file;
        private VCFHeader header = null;
        BcfToolsViewImpl(final String file) {
            this.file= file;
        }
        
        public String getFile() {
            return this.file;
        }
        public VCFHeader getHeader() {
            if(this.header==null) {
                try(VCFIterator r=open()) {
                    //should load the vcf header
                    }
                }
            return this.header;
            }
        public VCFIterator open() {
            return _query(null);
            }
        public VCFIterator query(final Locatable locatable) {
            return _query(locatable);
            }
        
        private VCFIterator _query(final Locatable locatable) {
            final List<String> cmd = getCommandViewPreamble();
            cmd.add(getFile());
            if(locatable!=null) {
                cmd.add(locatable.getContig()+":"+locatable.getStart()+":"+locatable.getEnd());
                }
            return _create(cmd);
            }
        

        }

    
    @Override
    public CloseableTribbleIterator<VariantContext> query(Locatable locus) throws IOException {
        return FeatureReader.super.query(locus);
        } 
    
    @Override
    public CloseableTribbleIterator<VariantContext> query(String chr, int start, int end)
            throws IOException {
        return query(new Interval(chr, start, end));
    }

    @Override
    public CloseableTribbleIterator<VariantContext> iterator() throws IOException {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void close() throws IOException {
       //do nothing
    }

    @Override
    public List<String> getSequenceNames() {
        // TODO Auto-generated method stub
        return null;
    }
   
    
    private class BcfToolsVcfIterator implements VCFIterator  {
        final Process proc;
        final VCFIterator delegate;
        public BcfToolsVcfIterator(VCFIterator delegate,final Process proc) {
            this.delegate = delegate;
            this.proc = proc;
            }
        @Override
        public void close() {
            this.delegate.close();
            this.proc.destroy();
            }
        }
    
    private static class Iter implements CloseableTribbleIterator<VariantContext> {
        @Override
        public void close() {
            
        }

    }
}
