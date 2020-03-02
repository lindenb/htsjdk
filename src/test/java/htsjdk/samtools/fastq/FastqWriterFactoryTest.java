package htsjdk.samtools.fastq;

import htsjdk.HtsjdkTest;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

public final class FastqWriterFactoryTest extends HtsjdkTest {

    @DataProvider(name="conditions")
    Object[][] getFastqs() {
        return new Object[][] {
            {false,false,false},
            {true,false,false},
            {false,false,true},
            {false,true,false},
        };
    }
    
    @Test(dataProvider="conditions")
    public void testCreateFastq(boolean compressed, boolean createMd5, boolean async) throws IOException {
        final Path out= Files.createTempFile("tmp.", ".fastq" +(compressed?".gz":""));
        final FastqRecord rec = new FastqRecord("FQ1", "ATGC", "+", "####");
        
        try(FastqWriter qw = new FastqWriterFactory().setUseAsyncIo(async).setCreateMd5(createMd5).newWriter(out)) {
            qw.write(rec);
        }
       
        try(FastqReader fq= new FastqReader(out)) {
            Assert.assertTrue(fq.hasNext());
            final FastqRecord rec2 = fq.next();
            Assert.assertEquals(rec, rec2);
            Assert.assertFalse(fq.hasNext());
        }
        
        Files.delete(out);
        if(createMd5)  Files.delete(out.resolve(".md5"));
    }
}
