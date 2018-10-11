package org.mskcc.netbox.test.genomic;

import junit.framework.TestCase;
import org.mskcc.netbox.genomic.ByteProfileData;
import org.mskcc.netbox.genomic.GeneticAlterationType;
import org.mskcc.netbox.genomic.util.CaseSetReader;
import org.mskcc.netbox.genomic.util.MutationReader;
import org.mskcc.netbox.genomic.util.TabDelimReader;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;

/**
 * JUnit Tests for Genomic Readers.
 */
public class TestReaders extends TestCase {

    /**
     * Tests the Genomic Readers.
     *
     * @throws IOException IO Error.
     */
    public final void testProfileData() throws IOException {
        CaseSetReader caseReader = new CaseSetReader(new File("testData/case_ids_test.txt"));
        HashSet<String> caseIdSet = caseReader.getCaseIdSet();

        TabDelimReader cnaReader = new TabDelimReader(new File("testData/cna_test.txt"),
                GeneticAlterationType.COPY_NUMBER_ALTERATION, caseIdSet, false);
        ByteProfileData byteProfileData = cnaReader.getByteProfile();
        byte byteValue = byteProfileData.getValue("9636", "TCGA-02-0001");
        assertEquals(1, byteValue);
        byteValue = byteProfileData.getValue("9636", "TCGA-06-0646");
        assertEquals(1, byteValue);
        assertEquals(9, byteProfileData.getGeneList().size());

        MutationReader mutationReader =
                new MutationReader(new File("testData/mutation_test.txt"), caseIdSet, false);
        byteProfileData = mutationReader.getByteProfile();
        assertEquals(0, byteProfileData.getValue("6656", "TCGA-02-0055"));
        assertEquals(1, byteProfileData.getValue("7490", "TCGA-02-0089"));
        assertEquals(1, byteProfileData.getValue("3371", "TCGA-06-0129"));
        assertEquals(95, byteProfileData.getGeneList().size());
        assertEquals(85, byteProfileData.getCaseIdList().size());
    }
}
