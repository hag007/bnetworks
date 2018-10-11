package org.mskcc.netbox.genomic;

import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.ProgressMonitor;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

/**
 * Utility Class for Summarizing Profile Data.
 */
public final class ProfileDataSummary {
    private HashMap<String, Double> geneAlteredMap = new HashMap<String, Double>();
    private ArrayList<GeneWithScore> geneAlteredList = new ArrayList<GeneWithScore>();
    private HashMap<String, Boolean> caseAlteredMap = new HashMap<String, Boolean>();
    private ArrayList<String> observedCaseList;
    private ArrayList<String> observedGeneList;
    private ByteProfileData profileData;
    private ProgressMonitor pMonitor;

    /**
     * Constructor.
     *
     * @param data Profile Data Object.
     */
    public ProfileDataSummary(ByteProfileData data) {
        this.pMonitor = ProgressMonitor.getInstance();
        this.profileData = data;
        observedGeneList = data.getGeneList();
        observedCaseList = data.getCaseIdList();
        geneAlteredList = determineFrequencyOfGeneAlteration();
    }

    /**
     * Gets the embedded ByteProfileData Data Object.
     *
     * @return ByteProfileData Object.
     */
    public ByteProfileData getProfileData() {
        return this.profileData;
    }

    /**
     * Gets the List of Observed Genes.
     *
     * @return List of Observed Genes.
     */
    public ArrayList<String> getObservedGeneList() {
        return observedGeneList;
    }

    /**
     * Gets the List of Observed Cases.
     *
     * @return List of Observed Cases.
     */
    public ArrayList<String> getObservedCaseList() {
        return observedCaseList;
    }

    /**
     * Gets the list of gene ordered by frequency of alteration.
     *
     * @return ArrayList of GeneWithScore Objects.
     */
    public ArrayList<GeneWithScore> getGeneFrequencyList() {
        return geneAlteredList;
    }

    /**
     * Gets percent of cases where gene X is altered.
     *
     * @param gene gene symbol.
     * @return percentage value.
     */
    public double getPercentCasesWhereGeneIsAltered(String gene) {
        return geneAlteredMap.get(gene);
    }

    /**
     * Gets percent of cases where gene set X is altered.
     *
     * @param gList Gene List.
     * @return percentage value.
     */
    public double getPercentCasesWhereGeneSetisAltered(ArrayList<String> gList) {
        int numCasesAffected = 0;

        //  Iterate through all cases
        for (String caseId : observedCaseList) {
            boolean caseIsAltered = false;

            //  Iterate through all gene lists
            for (String gene : gList) {
                boolean isAltered = isGeneAltered(gene, caseId);
                if (isAltered) {
                    caseIsAltered = true;
                }
            }
            if (caseIsAltered) {
                numCasesAffected++;
            }
        }
        return numCasesAffected / (double) observedCaseList.size();
    }

    /**
     * Gets percent of cases where gene set X is altered.
     *
     * @param gList Gene List.
     * @param caseSet Case Set to Test.
     * @return percentage value.
     */
    public double getPercentCasesWhereGeneSetisAltered(ArrayList<String> gList,
                                                       ArrayList<String> caseSet) {
        int numCasesAffected = 0;

        //  Iterate through all cases
        for (String caseId : caseSet) {
            boolean caseIsAltered = false;

            //  Iterate through all gene lists
            for (String gene : gList) {
                boolean isAltered = isGeneAltered(gene, caseId);
                if (isAltered) {
                    caseIsAltered = true;
                }
            }
            if (caseIsAltered) {
                numCasesAffected++;
            }
        }
        return numCasesAffected / (double) caseSet.size();
    }

    /**
     * Does this case contain an aleration in at least one gene in the set.
     *
     * @param caseId case ID.
     * @return true or false.
     */
    public boolean isCaseAltered(String caseId) {
        return caseAlteredMap.get(caseId);
    }

    /**
     * Determines if the gene X is altered in case Y.
     *
     * @param gene   gene symbol.
     * @param caseId case Id.
     * @return true or false.
     */
    public boolean isGeneAltered(String gene, String caseId) {
        byte value = profileData.getValue(gene, caseId);
        if (value == 0) {
            return false;
        } else {
            return true;
        }
    }


    /**
     * Determines frequency with which each gene is altered.
     */
    private ArrayList<GeneWithScore> determineFrequencyOfGeneAlteration() {
        ArrayList<GeneWithScore> localGeneList = new ArrayList<GeneWithScore>();
        //  First, determine frequency with which each gene is altered.
        //  Iterate through all genes.
        pMonitor.setMaxValue(observedGeneList.size());
        pMonitor.setCurrentMessage("Determining frequency of alteration for each observed gene.");
        for (String gene : observedGeneList) {
            int numSamplesWhereGeneIsAltered = 0;

            //  Iterate through all cases.
            for (String caseId : observedCaseList) {

                //  Determine if gene is altered in this case
                boolean isAltered = isGeneAltered(gene, caseId);

                //  If gene is altered in this case, increment counter.
                if (isAltered) {
                    numSamplesWhereGeneIsAltered++;
                }
            }
            CommandLineUtil.showProgress(pMonitor);
            pMonitor.incrementCurValue();
            double percent = numSamplesWhereGeneIsAltered / (double) observedCaseList.size();
            geneAlteredMap.put(gene, percent);
            GeneWithScore geneWithScore = new GeneWithScore();
            geneWithScore.setGene(gene);
            geneWithScore.setScore(percent);
            localGeneList.add(geneWithScore);
        }

        //  Sort genes, based on frequency of alteration.
        Collections.sort(localGeneList, new GeneWithScoreComparator());
        return localGeneList;
    }
}

/**
 * Ranks Genes based on core.
 */
class GeneWithScoreComparator implements Comparator {

    /**
     * Ranks genes based on score.
     *
     * @param o0 GeneWithScore0
     * @param o1 GeneWithScore1
     * @return integer value.
     */
    public int compare(Object o0, Object o1) {
        GeneWithScore gene0 = (GeneWithScore) o0;
        GeneWithScore gene1 = (GeneWithScore) o1;
        double score0 = gene0.getScore();
        double score1 = gene1.getScore();
        if (score0 == score1) {
            return gene0.getGene().compareTo(gene1.getGene());
        } else {
            if (score0 < score1) {
                return +1;
            } else {
                return -1;
            }
        }
    }
}
