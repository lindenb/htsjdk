package htsjdk.tribble.gene;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import htsjdk.tribble.annotation.Strand;

public abstract class AbstractGeneBaseData {
    private final String contig;
    private final String source;
    private final String type;
    private final int start;
    private final int end;
    private final double score;
    private final Strand strand;
    private final int phase;
    private final Map<String, List<String>> attributes;

    protected AbstractGeneBaseData(final String contig, final String source, final String type,
                 final int start, final int end, final Double score, final Strand strand, final int phase,
                 final Map<String, List<String>> attributes) {
        this.contig = contig;
        this.source = source;
        this.type = type;
        this.start = start;
        this.end = end;
        this.score = score;
        this.phase = phase;
        this.strand = strand;
        this.attributes = copyAttributesSafely(attributes);
    	}

    private static Map<String, List<String>> copyAttributesSafely(final Map<String, List<String>> attributes) {
        final Map<String, List<String>> modifiableDeepMap = new LinkedHashMap<>();

        for (final Map.Entry<String, List<String>> entry : attributes.entrySet()) {
            final List<String> unmodifiableDeepList = Collections.unmodifiableList(new ArrayList<>(entry.getValue()));
            modifiableDeepMap.put(entry.getKey(), unmodifiableDeepList);
        }

        return Collections.unmodifiableMap(modifiableDeepMap);
    }

    @Override
    public abstract boolean equals(Object other);

    @Override
    public abstract int hashCode();

    public String getContig() {
        return contig;
    }

    public String getSource() {
        return source;
    }

    public String getType() {
        return type;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public double getScore() {
        return score;
    }

    public Strand getStrand() {
        return strand;
    }

    public int getPhase() {
        return phase;
    }

    public Map<String, List<String>> getAttributes() {
        return attributes;
    }

    public List<String> getAttribute(final String key) {
        return attributes.getOrDefault(key, Collections.emptyList());
    }

}
