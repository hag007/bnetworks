package org.mskcc.netbox.model;

import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.NamedQueries;
import javax.persistence.NamedQuery;

/**
 * Interaction Object.
 */
@Entity
@NamedQueries({
        @NamedQuery(name = "org.mskcc.netbox.deleteAllInteractions",
                query = "delete from Interaction"),
        @NamedQuery(name = "org.mskcc.netbox.getAllInteractions",
                query = "from Interaction as interaction"),
        @NamedQuery(name = "org.mskcc.netbox.getInteractionsByGeneSymbol",
                query = "from Interaction as interaction where interaction.geneA "
                        + "= :geneSymbol OR interaction.geneB = :geneSymbol")
})
public final class Interaction {
    private String geneA;
    private String geneB;
    private String interactionType;
    private String experimentTypes;
    private String pmids;
    private String source;

    /**
     * The synthetic database key associated with this interaction once it is
     * persisted to a database.
     */
    @Id
    @GeneratedValue(strategy = GenerationType.AUTO)
    private Integer id;

    /**
     * Gets symbol for Gene A.
     *
     * @return symbol for Gene A.
     */
    public String getGeneA() {
        return geneA;
    }

    /**
     * Sets symbol for Gene A.
     *
     * @param a symbol for Gene A.
     */
    public void setGeneA(String a) {
        this.geneA = a;
    }

    /**
     * Gets symbol for Gene B.
     *
     * @return symbol for Gene B.
     */
    public String getGeneB() {
        return geneB;
    }

    /**
     * Sets symbol for Gene B.
     *
     * @param b symbol for Gene B.
     */
    public void setGeneB(String b) {
        this.geneB = b;
    }

    /**
     * Gets the Interaction Type.
     *
     * @return interaction type.
     */
    public String getInteractionType() {
        return interactionType;
    }

    /**
     * Sets the Interaction Type.
     *
     * @param type interaction type.
     */
    public void setInteractionType(String type) {
        this.interactionType = type;
    }

    /**
     * Gets the Experiment Types.
     *
     * @return experiment types.
     */
    public String getExperimentTypes() {
        return experimentTypes;
    }

    /**
     * Sets the Experiment Types.
     *
     * @param expTypes experiment types.
     */
    public void setExperimentTypes(String expTypes) {
        this.experimentTypes = expTypes;
    }

    /**
     * Gets the PMIDs.
     *
     * @return PMIDs.
     */
    public String getPmids() {
        return pmids;
    }

    /**
     * Sets the PMIDs.
     *
     * @param p PMIDs.
     */
    public void setPmids(String p) {
        this.pmids = p;
    }

    /**
     * Gets the data source, where this interaction comes from, e.g. REACTOME.
     *
     * @return data source.
     */
    public String getSource() {
        return source;
    }

    /**
     * Sets the data source, where this interaction comes from, e.g. REACTOME.
     *
     * @param s data source
     */
    public void setSource(String s) {
        this.source = s;
    }

    @Override
    /**
     * Overrides toString()
     */
    public String toString() {
        return "Interaction:  " + geneA + " " + interactionType + " " + geneB + ", " + source;
    }

    /**
     * Provides a Cytoscape SIF Version of this Interaction.
     *
     * @return SIF Text.
     */
    public String toSif() {
        return geneA + " " + interactionType + " " + geneB;
    }
}
