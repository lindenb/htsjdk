/**
 *
 *
 */
package htsjdk.tribble.bgen;
import htsjdk.samtools.util.CloseableIterator;
import java.io.IOException;
import java.nio.file.Paths;
import htsjdk.tribble.index.*;
import htsjdk.tribble.index.linear.*;
import java.nio.file.Path;

public class BGenTest {

private void writeIndex( Path path) throws Exception {
	LinearIndex index= IndexFactory.createLinearIndex(path,new BGenCodec());
	index.write(Paths.get(path.toString()+".idx"));
	}

private void read1(String filename) throws Exception {
System.err.println(filename);
	final Path path = Paths.get(filename);
	try(BGenFileReader r= new BGenFileReader(path)) {
		try(CloseableIterator<BGenFeature> iter=r.iterator()) {
			while(iter.hasNext()) {
			System.err.println(filename);
				iter.next();
				}
			}
		}
	writeIndex(path);
	}

private void doWork() {
	try {
		//read1("/home/lindenb/src/htsjdk/src/test/resources/htsjdk/tribble/bgen/complex.bgen");
		read1("/home/lindenb/src/htsjdk/src/test/resources/htsjdk/tribble/bgen/haplotypes.bgen");
		}
	catch(Throwable err) {
		err.printStackTrace();
		}
	}

public static void main(String[] args ) {
	try {
		new BGenTest().doWork();
	} catch(Throwable err) {
	err.printStackTrace();
	}
	}
}

