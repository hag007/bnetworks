package org.mskcc.netbox.model;

import java.util.ArrayList;

/**
 * Encapsualtes Information Regarding Gene Sets.
 *
 * @author Ethan Cerami.
 */
public final class GeneSet {
    private String name;
    private ArrayList<Gene> geneList;

    /**
     * Gets Gene Set Name.
     * @return gene set name.
     */
    public String getName() {
        return name;
    }

    /**
     * Sets Gene Set Name.
     * @param geneSetName gene set name.
     */
    public void setName(String geneSetName) {
        this.name = geneSetName;
    }

    /**
     * Gets the Gene List.
     * @return gene list.
     */
    public ArrayList<Gene> getGeneList() {
        return geneList;
    }

    /**
     * Sets the Gene List.
     * @param list gene list.
     */
    public void setGeneList(ArrayList<Gene> list) {
        this.geneList = list;
    }
}
