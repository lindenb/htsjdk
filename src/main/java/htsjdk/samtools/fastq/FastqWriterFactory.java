package htsjdk.samtools.fastq;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.AbstractAsyncWriter;
import java.io.File;
import java.nio.file.Path;

/**
 * Factory class for creating FastqWriter objects.
 *
 * @author Tim Fennell
 */
public class FastqWriterFactory {
    private boolean useAsyncIo = Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS;
    private boolean createMd5  = Defaults.CREATE_MD5;

    /** Sets whether or not to use async io (i.e. a dedicated thread per writer. */
    public FastqWriterFactory setUseAsyncIo(final boolean useAsyncIo) { this.useAsyncIo = useAsyncIo; return this;}

    /** @return is async io used */
    public boolean isUseAsyncIo() { return this.useAsyncIo;}
    
    /** If true, compute MD5 and write appropriately-named file when file is closed. */
    public FastqWriterFactory setCreateMd5(final boolean createMd5) { this.createMd5 = createMd5; return this;}

    /** @return is a MD5 checksum computed */
    public boolean isCreateMd5() { return this.createMd5;}
    
    public FastqWriter newWriter(final File out) {
        return newWriter(out.toPath());
    }
    
    public FastqWriter newWriter(final Path out) {
        final FastqWriter writer = new BasicFastqWriter(out, createMd5);
        if (useAsyncIo) {
            return new AsyncFastqWriter(writer, AbstractAsyncWriter.DEFAULT_QUEUE_SIZE);
        }
        else {
            return writer;
        }
    }
}
