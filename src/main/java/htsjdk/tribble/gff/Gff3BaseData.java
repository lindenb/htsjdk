package htsjdk.tribble.gff;

import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gene.AbstractGeneBaseData;

import java.util.Collections;
import java.util.List;
import java.util.Map;

public class Gff3BaseData extends AbstractGeneBaseData {
    private final String id;
    private final String name;
    private final List<String> aliases;
    private final int hashCode;

    public Gff3BaseData(final String contig, final String source, final String type,
                 final int start, final int end, final Double score, final Strand strand, final int phase,
                 final Map<String, List<String>> attributes) {
    	super(contig,source,type,start,end,score,strand,phase,attributes);
        this.id = Gff3Codec.extractSingleAttribute(attributes.get(Gff3Constants.ID_ATTRIBUTE_KEY));
        this.name = Gff3Codec.extractSingleAttribute(attributes.get(Gff3Constants.NAME_ATTRIBUTE_KEY));
        this.aliases = attributes.getOrDefault(Gff3Constants.ALIAS_ATTRIBUTE_KEY, Collections.emptyList());
        this.hashCode = computeHashCode();
    }


    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if(!other.getClass().equals(Gff3BaseData.class)) {
            return false;
        }

        final Gff3BaseData otherBaseData = (Gff3BaseData) other;
        boolean ret = otherBaseData.getContig().equals(getContig()) &&
                otherBaseData.getSource().equals(getSource()) &&
                otherBaseData.getType().equals(getType()) &&
                otherBaseData.getStart() == getStart() &&
                otherBaseData.getEnd() == getEnd() &&
                ((Double)otherBaseData.getScore()).equals(this.getScore()) &&
                otherBaseData.getPhase() == getPhase() &&
                otherBaseData.getStrand().equals(getStrand()) &&
                otherBaseData.getAttributes().equals(getAttributes());
        if (getId() == null) {
            ret = ret && otherBaseData.getId() == null;
        } else {
            ret = ret && otherBaseData.getId() != null && otherBaseData.getId().equals(getId());
        }

        if (getName() == null) {
            ret = ret && otherBaseData.getName() == null;
        } else {
            ret = ret && otherBaseData.getName() != null && otherBaseData.getName().equals(getName());
        }

        ret = ret && otherBaseData.getAliases().equals(getAliases());

        return ret;
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    private int computeHashCode() {
        int hash = getContig().hashCode();
        hash = 31 * hash + getSource().hashCode();
        hash = 31 * hash + getType().hashCode();
        hash = 31 * hash + getStart();
        hash = 31 * hash + getEnd();
        hash = 31 * hash + Double.hashCode(getScore());
        hash = 31 * hash + getPhase();
        hash = 31 * hash + getStrand().hashCode();
        hash = 31 * hash + getAttributes().hashCode();
        if (getId() != null) {
            hash = 31 * hash + getId().hashCode();
        }

        if (getName() != null) {
            hash = 31 * hash + getName().hashCode();
        }

        hash = 31 * hash + aliases.hashCode();

        return hash;
    }

    public String getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public List<String> getAliases() {
        return aliases;
    }
}
