/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package htsjdk.variant.variantcontext;


import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.vcf.VCFConstants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * This class encompasses all the basic information about a genotype.  It is immutable.
 *
 * @author Mark DePristo
 */
abstract class AbstractGenotype implements Genotype {
    public static final long serialVersionUID = 1L;

    private final String sampleName;
    private GenotypeType type = null;
    private final String filters;

    protected AbstractGenotype(final String sampleName, final String filters) {
        this.sampleName = sampleName;
        this.filters = filters == null || filters.isEmpty() ? null : filters;
    }

    /**
     * Returns the name associated with this sample.
     *
     * @return a non-null String
     */
    @Override
    public String getSampleName() {
        return sampleName;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // The type of this genotype
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * @return the high-level type of this sample's genotype
     */
    @Override
    public GenotypeType getType() {
        if ( type == null ) {
            type = determineType();
        }
        return type;
    }

    /**
     * Internal code to determine the type of the genotype from the alleles vector
     * @return the type
     */
    protected GenotypeType determineType() {
        // TODO -- this code is slow and could be optimized for the diploid case
        final List<Allele> alleles = getAlleles();
        if ( alleles.isEmpty() )
            return GenotypeType.UNAVAILABLE;

        boolean sawNoCall = false, sawMultipleAlleles = false;
        Allele observedAllele = null;

        for ( final Allele allele : alleles ) {
            if ( allele.isNoCall() )
                sawNoCall = true;
            else if ( observedAllele == null )
                observedAllele = allele;
            else if ( !allele.equals(observedAllele) )
                sawMultipleAlleles = true;
        }

        if ( sawNoCall ) {
            if ( observedAllele == null )
                return GenotypeType.NO_CALL;
            return GenotypeType.MIXED;
        }

        if ( observedAllele == null )
            throw new IllegalStateException("BUG: there are no alleles present in this genotype but the alleles list is not null");

        return sawMultipleAlleles ? GenotypeType.HET : observedAllele.isReference() ? GenotypeType.HOM_REF : GenotypeType.HOM_VAR;
    }


    // ---------------------------------------------------------------------------------------------------------
    //
    // Many different string representations
    //
    // ---------------------------------------------------------------------------------------------------------


    /**
     * Return a VCF-like string representation for the alleles of this genotype.
     *
     * If ignoreRefState is true, will not append the reference * marker on the alleles.
     *
     * @return a string representing the genotypes, or null if the type is unavailable.
     */
    public String getGenotypeString(boolean ignoreRefState) {
        if ( getPloidy() == 0 )
            return "NA";

        // Notes:
        // 1. Make sure to use the appropriate separator depending on whether the genotype is phased
        final String separator = isPhased() ? PHASED_ALLELE_SEPARATOR : UNPHASED_ALLELE_SEPARATOR;
        // 2. If ignoreRefState is true, then we want just the bases of the Alleles (ignoring the '*' indicating a ref Allele)
        if (ignoreRefState) {
          return ParsingUtils.join(separator, getAlleleStrings());
        }
        // 3. So that everything is deterministic with regards to integration tests, we sort Alleles (when the genotype isn't phased, of course)
        List<Allele> alleles = isPhased() ? getAlleles() : ParsingUtils.sortList(getAlleles());
        return ParsingUtils.join(separator, alleles);
    }

    /**
     * Utility that returns a list of allele strings corresponding to the alleles in this sample
     * @return
     */
    protected List<String> getAlleleStrings() {
        final List<String> al = new ArrayList<String>(getPloidy());
        for ( Allele a : getAlleles() )
            al.add(a.getBaseString());

        return al;
    }

    @Override
    public String toString() {
        return String.format("[%s %s%s%s%s%s%s%s]",
                getSampleName(),
                getGenotypeString(false),
                toStringIfExists(VCFConstants.GENOTYPE_QUALITY_KEY, getGQ()),
                toStringIfExists(VCFConstants.DEPTH_KEY, getDP()),
                toStringIfExists(VCFConstants.GENOTYPE_ALLELE_DEPTHS, getAD()),
                toStringIfExists(VCFConstants.GENOTYPE_PL_KEY, getPL()),
                toStringIfExists(VCFConstants.GENOTYPE_FILTER_KEY, getFilters()),
                sortedString(getExtendedAttributes()));
    }

    @Override
    public String toBriefString() {
        return String.format("%s:Q%d", getGenotypeString(false), getGQ());
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // get routines for extended attributes
    //
    // ---------------------------------------------------------------------------------------------------------

    /**
     * Returns the filter string associated with this Genotype.
     *
     * @return If this result == null, then the genotype is considered PASSing filters
     *   If the result != null, then the genotype has failed filtering for the reason(s)
     *   specified in result.  To be reference compliant multiple filter field
     *   string values can be encoded with a ; separator.
     */
    @Override
    public final String getFilters() {
        return filters;
    }


    // ------------------------------------------------------------------------------
    //
    // private utilities
    //
    // ------------------------------------------------------------------------------

    /**
     * a utility method for generating sorted strings from a map key set.
     * @param c the map
     * @param <T> the key type
     * @param <V> the value type
     * @return a sting, enclosed in {}, with comma seperated key value pairs in order of the keys
     */
    protected static <T extends Comparable<T>, V> String sortedString(Map<T, V> c) {

        // NOTE -- THIS IS COPIED FROM GATK UTILS TO ALLOW US TO KEEP A SEPARATION BETWEEN THE GATK AND VCF CODECS
        final List<T> t = new ArrayList<T>(c.keySet());
        Collections.sort(t);

        final List<String> pairs = new ArrayList<String>();
        for (final T k : t) {
            pairs.add(k + "=" + c.get(k));
        }

        return pairs.isEmpty() ? "" : " {" + ParsingUtils.join(", ", pairs.toArray(new String[pairs.size()])) + "}";
    }

    /**
     * Returns a display name for field name with value v if this isn't -1.  Otherwise returns ""
     * @param name of the field ("AD")
     * @param v the value of the field, or -1 if missing
     * @return a non-null string for display if the field is not missing
     */
    protected final static String toStringIfExists(final String name, final int v) {
        return v == -1 ? "" : " " + name + " " + v;
    }

    /**
     * Returns a display name for field name with String value v if this isn't null.  Otherwise returns ""
     * @param name of the field ("FT")
     * @param v the value of the field, or null if missing
     * @return a non-null string for display if the field is not missing
     */
    protected final static String toStringIfExists(final String name, final String v) {
        return v == null ? "" : " " + name + " " + v;
    }

    /**
     * Returns a display name for field name with values vs if this isn't null.  Otherwise returns ""
     * @param name of the field ("AD")
     * @param vs the value of the field, or null if missing
     * @return a non-null string for display if the field is not missing
     */
    protected final static String toStringIfExists(final String name, final int[] vs) {
        if ( vs == null )
            return "";
        else {
            StringBuilder b = new StringBuilder();
            b.append(' ').append(name).append(' ');
            for ( int i = 0; i < vs.length; i++ ) {
                if ( i != 0 ) b.append(',');
                b.append(vs[i]);
            }
            return b.toString();
        }
    }

    /**
     * Does the attribute map have a mapping involving a forbidden key (i.e.,
     * one that's managed inline by this Genotypes object?
     *
     * @param attributes the extended attributes key
     * @return
     */
    protected final static boolean hasForbiddenKey(final Map<String, Object> attributes) {
        for ( final String forbidden : PRIMARY_KEYS)
            if ( attributes.containsKey(forbidden) )
                return true;
        return false;
    }

    protected final static boolean isForbiddenKey(final String key) {
        return PRIMARY_KEYS.contains(key);
    }
}
