package org.mskcc.netbox.genomic;

import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.GlobalConfig;
import org.mskcc.netbox.util.ProgressMonitor;

import java.util.ArrayList;

/**
 * Merges Multiple ByteProfiles Into One.
 *
 * @author Ethan Cerami.
 */
public final class ByteProfileMerger {
    private ByteProfileData mergedProfile;
    private ProgressMonitor pMonitor;

    /**
     * Constructor.
     *
     * @param profileList ArrayList of ByteProfileData Objects.
     */
    public ByteProfileMerger(ArrayList<ByteProfileData> profileList) {
        pMonitor = ProgressMonitor.getInstance();
        pMonitor.setCurrentMessage("Merging Profile Data");

        //  Create Union of all Cases and all Genes
        ArrayList<String> caseList = new ArrayList<String>();
        ArrayList<String> geneList = new ArrayList<String>();
        createUnion(profileList, caseList, geneList);

        //  Perform the actual merge
        mergedProfile = new ByteProfileData(GeneticAlterationType.MERGED_PROFILE,
                geneList, caseList);
        mergeProfiles(profileList, caseList, geneList);
    }

    /**
     * Gets the new merged profile data object.
     *
     * @return ByteProfileData Object.
     */
    public ByteProfileData getMergedProfile() {
        return mergedProfile;
    }

    /**
     * Perform the merge.
     */
    private void mergeProfiles(ArrayList<ByteProfileData> profileList,
                               ArrayList<String> caseList, ArrayList<String> geneList) {
        //  Iterate through all genes
        pMonitor.setMaxValue(geneList.size());

        for (String gene : geneList) {

            //  Iterate through all cases
            for (String caseId : caseList) {

                //  Determine status of gene X in caseId Y.
                byte status = isAltered(profileList, gene, caseId);

                mergedProfile.setValue(gene, caseId, status);
            }
            CommandLineUtil.showProgress(pMonitor);
            pMonitor.incrementCurValue();
        }
    }

    /**
     * Determines the alteration status of gene X in case Y.
     */
    private byte isAltered(ArrayList<ByteProfileData> profileList,
                           String gene, String caseId) {

        //  Iterate through all profiles
        byte isAltered = 0;
        for (ByteProfileData data : profileList) {

            //  Get the data value and alteration type for gene X in caseId Y in profile Z
            byte value = data.getValue(gene, caseId);
            GeneticAlterationType alterationType = data.getAlterationType();

            //  Handle Copy Number Changes
            if (alterationType.equals(GeneticAlterationType.COPY_NUMBER_ALTERATION)) {
                if (GlobalConfig.getInstance().includeLowLevelCnaChanges()) {
                    if (value != 0) {
                        isAltered = 1;
                    }
                } else {
                    if (value == -2 || value == 2) {
                        isAltered = 1;
                    }
                }
            } else if (alterationType.equals(GeneticAlterationType.MUTATION)) {
                if (value == 1) {
                    isAltered = 1;
                }
            } else if (alterationType.equals(GeneticAlterationType.MRNA_EXPRESSION)) {
                if (value == 1) {
                    isAltered = 1;
                }
            }
        }
        return isAltered;
    }

    /**
     * Creates the Union of all Cases and the Union of all Genes.
     */
    private void createUnion(ArrayList<ByteProfileData> profileList,
                             ArrayList<String> caseIdList, ArrayList<String> geneList) {

        //  Iterate through all profiles
        for (ByteProfileData data : profileList) {

            //  Get the case list and the gene list
            ArrayList<String> currentCaseList = data.getCaseIdList();
            ArrayList<String> currentGeneList = data.getGeneList();

            //  Conditionally add each new case to the global case list
            for (String currentCaseId : currentCaseList) {
                if (!caseIdList.contains(currentCaseId)) {
                    caseIdList.add(currentCaseId);
                }
            }

            //  Conditionally add each new gene to the global gene list
            for (String currentGene : currentGeneList) {
                if (!geneList.contains(currentGene)) {
                    geneList.add(currentGene);
                }
            }
        }
    }
}
