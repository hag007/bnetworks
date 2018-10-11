package org.mskcc.netbox.genomic;

import org.mskcc.netbox.util.Formatter;

/**
 * Encapsulates a Gene with a Score.
 * Score is used to rank genes.
 */
public final class GeneWithScore {
    private String gene;
    private int entrezGeneId;
    private double score;

    /**
     * Gets the Gene Symbol.
     *
     * @return gene symbol.
     */
    public String getGene() {
        return gene;
    }

    /**
     * Sets the Gene Symbol.
     *
     * @param g gene symbol.
     */
    public void setGene(String g) {
        gene = g;
    }

    /**
     * Gets Entrez Gene Id.
     *
     * @return Entrez Gene Id.
     */
    public int getEntrezGeneId() {
        return entrezGeneId;
    }

    /**
     * Sets Entrez Gene Id.
     *
     * @param id Entrez Gene Id;
     */
    public void setEntrezGeneId(int id) {
        entrezGeneId = id;
    }

    /**
     * Gets the score.
     *
     * @return gene score.
     */
    public double getScore() {
        return score;
    }

    /**
     * Sets the gene score.
     *
     * @param s gene score.
     */
    public void setScore(double s) {
        score = s;
    }

    /**
     * Overrides toString().
     *
     * @return gene info.
     */
    public String toString() {
        return gene + " (Altered at frequency:  "
                + Formatter.getDecimalFormat().format(score) + ")";
    }
}
