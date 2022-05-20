package htsjdk.tribble.gff;


/**
 * Gff3 format spec is defined at https://github.com/The-Sequence-Ontology/Specifications/blob/31f62ad469b31769b43af42e0903448db1826925/gff3.md
 * Discontinuous features which are split between multiple lines in the gff files are implemented as separate features linked as "co-features"
 */
public interface Gff3Feature extends GeneFeature<Gff3Feature> {

}
