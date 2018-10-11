package org.mskcc.netbox.genomic.util;

import org.mskcc.netbox.genomic.ByteProfileData;
import org.mskcc.netbox.genomic.GeneticAlterationType;
import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.FileUtil;
import org.mskcc.netbox.util.ProgressMonitor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Reads Mutation Data into a ByteProfile Data Object.
 *
 * @author Ethan Cerami
 */
public final class MutationReader {
    private ProgressMonitor pMonitor;
    private ByteProfileData byteProfile;

    /**
     * Constructor.
     *
     * @param file          Mutation File.
     * @param caseIdSet     Case ID Set.
     * @param validateGenes Flag to Validate Genes.
     * @throws IOException IO Error.
     */
    public MutationReader(File file, HashSet<String> caseIdSet, boolean validateGenes)
            throws IOException {
        this.pMonitor = ProgressMonitor.getInstance();
        execute(file, caseIdSet, validateGenes);
    }

    /**
     * Gets the Byte Profile Data.
     *
     * @return ByteProfile Data.
     */
    public ByteProfileData getByteProfile() {
        return byteProfile;
    }

    private void execute(File file, HashSet<String> caseIdSet, boolean validateGenes)
            throws IOException {
        int numLines = FileUtil.getNumLines(file);

        //  Pass #1:  Determine Number of Genes
        //  We assume first line contains the header info and ignore it
        pMonitor.setCurrentMessage("Reading data from:  " + file.getAbsolutePath()
                + ", Pass #1");
        pMonitor.setMaxValue(numLines);
        FileReader reader = new FileReader(file);
        BufferedReader buf = new BufferedReader(reader);
        buf.readLine();

        HashSet<String> geneSet = new HashSet<String>();

        String line = buf.readLine();
        while (line != null) {
            if (line.trim().length() > 0) {
                String[] parts = line.split("\t");
                String entrezGeneId = parts[0];
                String id = ReaderUtil.getGeneId(entrezGeneId, validateGenes, pMonitor);
                if (id != null) {
                    geneSet.add(id);
                }
            }
            line = buf.readLine();
            pMonitor.incrementCurValue();
            CommandLineUtil.showProgress(pMonitor);
        }

        //  At this point, we can create the ByteProfile Data Object
        ArrayList<String> geneList = new ArrayList<String>();
        geneList.addAll(geneSet);
        ArrayList<String> caseIdList = new ArrayList<String>();
        caseIdList.addAll(caseIdSet);
        byteProfile = new ByteProfileData(GeneticAlterationType.MUTATION, geneList, caseIdList);

        //  Pass #2:  Read in the Actual Data
        pMonitor.setCurrentMessage("Reading data from:  " + file.getAbsolutePath()
                + ", Pass #2");
        pMonitor.setMaxValue(numLines);

        reader = new FileReader(file);
        buf = new BufferedReader(reader);

        buf.readLine();
        line = buf.readLine();
        while (line != null) {
            if (line.trim().length() > 0) {
                String[] parts = line.split("\t");
                String entrezGeneId = parts[0];
                String caseId = parts[1];
                String mutationCall = parts[2];
                byte mutCall = GenomicCall.getMutationCall(mutationCall);
                String id = ReaderUtil.getGeneId(entrezGeneId, validateGenes, pMonitor);

                //  Only include values that are from the selected case set.
                if (id != null && caseIdSet.contains(caseId)) {
                    byteProfile.setValue(id, caseId, mutCall);
                }
            }
            line = buf.readLine();
            pMonitor.incrementCurValue();
            CommandLineUtil.showProgress(pMonitor);
        }
    }
}
