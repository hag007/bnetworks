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
 * Reads Tab Delimited Data into a ByteProfile Data Object.
 *
 * @author Ethan Cerami
 */
public final class TabDelimReader {
    private ProgressMonitor pMonitor;
    private ByteProfileData byteProfile;
    private GeneticAlterationType alterationType;
    private HashSet<String> caseIdSet;
    private File file;
    private boolean validateGenes;
    private int numLines;
    private int startDataColumn;

    /**
     * Constructor.
     *
     * @param f                 Input File.
     * @param type              Genetic Alteration Type.
     * @param caseSet           HashSet of Case IDs to Include.
     * @param validateGenesFlag Flag to Validate Genes.
     * @throws IOException IO Error.
     */
    public TabDelimReader(File f, GeneticAlterationType type, HashSet<String> caseSet,
                          boolean validateGenesFlag) throws IOException {
        this.file = f;
        this.alterationType = type;
        this.caseIdSet = caseSet;
        this.validateGenes = validateGenesFlag;
        this.pMonitor = ProgressMonitor.getInstance();

        if (alterationType.equals(GeneticAlterationType.COPY_NUMBER_ALTERATION)) {
            startDataColumn = 2;
        } else {
            startDataColumn = 1;
        }

        execute();
    }

    private void execute() throws IOException {
        numLines = FileUtil.getNumLines(file);

        ArrayList<String> geneList = getGeneList();
        ArrayList<String> caseIdList = new ArrayList<String>();
        caseIdList.addAll(caseIdSet);

        byteProfile = new ByteProfileData(alterationType, geneList, caseIdList);

        pMonitor.setCurrentMessage("Reading data from:  " + file.getAbsolutePath()
                + ", Pass #2");

        pMonitor.setMaxValue(numLines);
        FileReader reader = new FileReader(file);
        BufferedReader buf = new BufferedReader(reader);
        String line = buf.readLine();

        //  We assume first line contains the header info with case IDs.
        String[] caseIds = line.split("\t");

        line = buf.readLine();
        while (line != null) {
            if (!line.startsWith("#")) {
                String[] parts = line.split("\t");
                String entrezGeneId = getEntrezGeneId(parts);
                String id = ReaderUtil.getGeneId(entrezGeneId, validateGenes, pMonitor);
                if (id != null) {
                    addRow(id, caseIds, parts);
                }
            }
            line = buf.readLine();
            pMonitor.incrementCurValue();
            CommandLineUtil.showProgress(pMonitor);
        }
    }

    /**
     * Gets the ByteProfile Data Object.
     *
     * @return ByteProfile Data Object.
     */
    public ByteProfileData getByteProfile() {
        return byteProfile;
    }

    private String getEntrezGeneId(String[] parts) {
        if (alterationType.equals(GeneticAlterationType.COPY_NUMBER_ALTERATION)) {
            return parts[1];
        } else {
            return parts[0];
        }
    }

    /**
     * Conducts First Pass at Data File to Determine Gene List.
     */
    private ArrayList<String> getGeneList()
            throws IOException {
        pMonitor.setCurrentMessage("Reading data from:  " + file.getAbsolutePath()
                + ", Pass #1");
        pMonitor.setMaxValue(numLines);

        ArrayList<String> geneList = new ArrayList<String>();
        FileReader reader = new FileReader(file);
        BufferedReader buf = new BufferedReader(reader);
        buf.readLine();
        String line = buf.readLine();
        while (line != null) {
            if (!line.startsWith("#")) {
                String[] parts = line.split("\t");
                String entrezGeneId = getEntrezGeneId(parts);
                String id = ReaderUtil.getGeneId(entrezGeneId, validateGenes, pMonitor);
                if (id != null) {
                    geneList.add(id);
                }
            }
            line = buf.readLine();
            pMonitor.incrementCurValue();
            CommandLineUtil.showProgress(pMonitor);
        }
        return geneList;
    }

    private void addRow(String geneId, String[] caseIds, String[] parts) {
        for (int i = startDataColumn; i < parts.length; i++) {
            String value = parts[i];
            byte byteValue = Byte.MIN_VALUE;
            if (alterationType.equals(GeneticAlterationType.MUTATION)) {
                byteValue = GenomicCall.getMutationCall(value);
            } else if (alterationType.equals(GeneticAlterationType.MRNA_EXPRESSION)) {
                byteValue = GenomicCall.getExpressionCall(value);
            } else {
                try {
                    byteValue = Byte.parseByte(value);
                } catch (NumberFormatException e) {
                    byteValue = Byte.MIN_VALUE;
                }
            }

            //  Only add values for select cases
            String caseId = caseIds[i];
            if (caseIdSet.contains(caseId)) {
                byteProfile.setValue(geneId, caseId, byteValue);
            }
        }
    }
}
