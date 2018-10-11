package org.mskcc.netbox.graph;

import org.mskcc.netbox.genomic.GeneWithScore;

import java.util.ArrayList;
import java.util.HashSet;

/**
 * Encapsulates Module Data.
 *
 * @author Ethan Cerami.
 */
public final class Module {
    private int moduleId;
    private ArrayList<GeneWithScore> geneList = new ArrayList<GeneWithScore>();
    private double pValueUnAdjusted;
    private double pValueFdrAdjusted;
    private double score;
    private String label;
    private HashSet<String> genesExplored = new HashSet<String>();

    /**
     * Gets the Module Identifier.
     *
     * @return module ID.
     */
    public int getModuleId() {
        return moduleId;
    }

    /**
     * Sets the Module Identifier.
     *
     * @param id module ID.
     */
    public void setModuleId(int id) {
        this.moduleId = id;
    }

    /**
     * Gets all Genes in the Module.
     *
     * @return List of GeneWithScore Objects.
     */
    public ArrayList<GeneWithScore> getGeneList() {
        return geneList;
    }

    /**
     * Sets all Genes in the Module.
     *
     * @param gList List of GeneWithScore Objects.
     */
    public void setGeneList(ArrayList<GeneWithScore> gList) {
        this.geneList = gList;
        if (geneList.size() > 0) {
            GeneWithScore gene = geneList.get(0);
            label = gene.getGene();
        }
    }

    /**
     * Gets the Unadjusted P-Value.
     *
     * @return unadjusted p-value.
     */
    public double getPValueUnAdjusted() {
        return pValueUnAdjusted;
    }

    /**
     * Sets the Unadjusted P-Value.
     *
     * @param p unadjusted p-value.
     */
    public void setPValueUnAdjusted(double p) {
        this.pValueUnAdjusted = p;
    }

    /**
     * Gets the FDR Adjuscted P-Value.
     *
     * @return FDR Adjusted P-Value.
     */
    public double getPValueFdrAdjusted() {
        return pValueFdrAdjusted;
    }

    /**
     * Sets the FDR Adjusted P-Value.
     *
     * @param p FDR Adjusted P-Value.
     */
    public void setPValueFdrAdjusted(double p) {
        this.pValueFdrAdjusted = p;
    }

    /**
     * Gets the Frequency of Alteration Score.
     *
     * @return score.
     */
    public double getScore() {
        return score;
    }

    /**
     * Sets the Frequency of Alteration Score.
     *
     * @param s score.
     */
    public void setScore(double s) {
        this.score = s;
    }

    /**
     * Gets the Module Label.
     *
     * @return Module Label.
     */
    public String getLabel() {
        return label;
    }

    /**
     * Gets the Genes explored during the greedy algorithm exploration.
     *
     * @return HashSet of Genes.
     */
    public HashSet<String> getGenesExplored() {
        return genesExplored;
    }

    /**
     * Sets the Genes explored during the greedy algorithm exploration.
     *
     * @param set HashSet of Genes.
     */
    public void setGenesExplored(HashSet<String> set) {
        this.genesExplored = set;
    }
}
