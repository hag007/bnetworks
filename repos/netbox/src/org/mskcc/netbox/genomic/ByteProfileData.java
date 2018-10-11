package org.mskcc.netbox.genomic;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Byte Profile Data Object.
 *
 * @author Ethan Cerami
 */
public final class ByteProfileData {
    private GeneticAlterationType alterationType;
    private byte[][] matrix;
    private ArrayList<String> geneList;
    private ArrayList<String> caseIdList;
    private HashMap<String, Integer> geneRowMap = new HashMap<String, Integer>();
    private HashMap<String, Integer> caseColMap = new HashMap<String, Integer>();


    /**
     * Constructor.
     *
     * @param type  Genetic Alteration Type.
     * @param gList List of Genes.
     * @param cList List of Case IDs.
     */
    public ByteProfileData(GeneticAlterationType type,
                           ArrayList<String> gList, ArrayList<String> cList) {
        this.alterationType = type;
        this.geneList = gList;
        this.caseIdList = cList;
        this.matrix = new byte[gList.size()][cList.size()];
        initRowCol(gList, cList);
    }

    /**
     * Constructor.
     *
     * @param type     Genetic Alteration Type.
     * @param gList    List of Genes.
     * @param caseList List of Case IDs.
     * @param m        2D Matrix of Byte Data.
     */
    public ByteProfileData(GeneticAlterationType type,
                           ArrayList<String> gList, ArrayList<String> caseList,
                           byte[][] m) {
        this.alterationType = type;
        this.geneList = gList;
        this.caseIdList = caseList;
        this.matrix = m;
        initRowCol(gList, caseList);
    }

    private void initRowCol(ArrayList<String> gList, ArrayList<String> caseList) {
        for (int i = 0; i < gList.size(); i++) {
            String gene = gList.get(i);
            geneRowMap.put(gene, i);
        }

        for (int i = 0; i < caseList.size(); i++) {
            String caseId = caseList.get(i);
            caseColMap.put(caseId, i);
        }
    }

    /**
     * Gets the Gene List.
     *
     * @return List of Gene Symbols.
     */
    public ArrayList<String> getGeneList() {
        return geneList;
    }

    /**
     * Gets the Case ID List.
     *
     * @return List of Case IDs.
     */
    public ArrayList<String> getCaseIdList() {
        return caseIdList;
    }

    /**
     * Gets the Value of the Specified Gene in the Specified Case ID.
     *
     * @param gene   Gene Symbol.
     * @param caseId Case ID.
     * @return byte value.
     */
    public byte getValue(String gene, String caseId) {
        try {
            int row = geneRowMap.get(gene);
            int col = caseColMap.get(caseId);
            return matrix[row][col];
        } catch (NullPointerException e) {
            return 0;
        }
    }

    /**
     * Sets the Values of the Specified Gene in the Specified Case ID.
     *
     * @param gene   Gene Symbol.
     * @param caseId Case ID.
     * @param b      byte value.
     */
    public void setValue(String gene, String caseId, byte b) {
        int row = geneRowMap.get(gene);
        int col = caseColMap.get(caseId);
        matrix[row][col] = b;
    }

    /**
     * Gets the Genetic Alteration Type.
     *
     * @return Genetic Alteration Type.
     */
    public GeneticAlterationType getAlterationType() {
        return alterationType;
    }
}
