package htsjdk.samtools.util;

/** MBean interface for ProgressLogger */
public interface ProgressLoggerMBean
    {
    /* the noun to use when logging, e.g. "Records, Variants, Loci" */
    public String getNoun();
    /* verb the verb to log, e.g. "Processed, Read, Written" */
    public String getVerb();
    /** Returns the count of records processed. */
    public long getCount();
    /** elapsed time */
    public String getElapsedTime();
    /** last record */
    public String getLastRecord();
    }

