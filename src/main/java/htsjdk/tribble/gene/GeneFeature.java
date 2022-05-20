package htsjdk.tribble.gene;

import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;

public interface GeneFeature<T extends GeneFeature> extends Feature {
    boolean isTopLevelFeature();


  default String getSource() {
      return getBaseData().getSource();
  }

  @Override
  default int getEnd() {
      return getBaseData().getEnd();
  }

  default Strand getStrand() {
      return getBaseData().getStrand();
  }

  default int getPhase() {
      return getBaseData().getPhase();
  }

  default String getType() {return getBaseData().getType();}

  @Override
  default String getContig() {
      return getBaseData().getContig();
  }

   @Override
  default int getStart() {
      return getBaseData().getStart();
  }


  default List<String> getAttribute(final String key) {
      return getBaseData().getAttribute(key);
  }

  default Map<String, List<String>> getAttributes() { return getBaseData().getAttributes();}

  default String getID() { return getBaseData().getId();}

  default String getName() { return getBaseData().getName();}

  default List<String> getAliases() { return getBaseData().getAliases();}

  default double getScore() { return getBaseData().getScore();}

  /**
   * Get BaseData object which contains all the basic information of the feature
   * @return
   */
  T getBaseData();

  /**
   * Gets set of parent features
   * @return set of parent features
   */
  Set<? extends T> getParents();

  /**
   * Gets set of features for which this feature is a parent
   * @return set of child features
   */
  Set<? extends T> getChildren();

  /**
   * Get set of all features this feature descends from, through chains of Parent attributes.  If Derives_From exists for this feature,
   * then only features along the inheritance path specified by the Derives_From attribute should be included as ancestors of this feature
   * @return set of ancestor features
   */
  Set<? extends T> getAncestors();

  /**
   * Get set of all features descended from this features, through chains of Parent attributes.  If Derives_From attribute exists for a feature,
   * it should only be included as a descendent of this feature if the inheritance path specified by its Derives_From attribute includes this feature
   * @return set of descendents
   */
  Set<? extends T> getDescendents();

  /**
   * Get set of co-features.  Co-features correspond to the other lines in the gff file that together make up a single discontinuous feature
   * @return set of co-features
   */
  Set<? extends T> getCoFeatures();

  /***
   * Flatten this feature and all descendents into a set of features.  The Derives_From attribute is respected if it exists
   * for this feature
   * @return set of this feature and all descendents
   */
  Set<? extends T> flatten();

  boolean hasParents();

  boolean hasChildren();

  boolean hasCoFeatures();

}
