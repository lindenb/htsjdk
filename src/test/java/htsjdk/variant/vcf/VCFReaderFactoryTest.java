/*
 * The MIT License
 *
 * Copyright (c) 2020 L'Institut-du-Thorax / Nantes
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package htsjdk.variant.vcf;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import htsjdk.HtsjdkTest;

public class VCFReaderFactoryTest  extends HtsjdkTest {
    @DataProvider(name = "pathsData")
    Object[][] pathsData() {
        final String TEST_DATA_DIR = "src/test/resources/htsjdk/variant/";
        return new Object[][]{
            {TEST_DATA_DIR + "test.vcf.bgz", true},
            {TEST_DATA_DIR + "test.vcf.bgz", false},
            {TEST_DATA_DIR + "VcfThatLacksAnIndex.vcf.gz", false},
            {TEST_DATA_DIR + "VcfThatLacksAnIndex.vcf",false},
            {TEST_DATA_DIR + "VcfThatLacksAnIndex but has a space.vcf",false}
            };
    }
    
    @Test(dataProvider = "pathsData")
    public void testOpenVCFReader(final String file, boolean withIndex) throws Exception {
        try(VCFReader r = VCFReaderFactory.newInstance().open(file)) {
            Assert.assertEquals(r.isQueryable(), withIndex);
        }
    }
}
