package org.mskcc.netbox.script;

import org.mskcc.netbox.genomic.ByteProfileData;
import org.mskcc.netbox.genomic.ByteProfileMerger;
import org.mskcc.netbox.genomic.GeneWithScore;
import org.mskcc.netbox.genomic.GeneticAlterationType;
import org.mskcc.netbox.genomic.ProfileDataSummary;
import org.mskcc.netbox.genomic.util.CaseSetReader;
import org.mskcc.netbox.genomic.util.MutationReader;
import org.mskcc.netbox.genomic.util.TabDelimReader;
import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.model.GeneSet;
import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.ProgressMonitor;
import org.mskcc.netbox.util.ReadGeneSets;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Analyzes Gene Sets.
 *
 * @author Ethan Cerami.
 */
public final class AnalyzeGeneSets {
    private static ProgressMonitor pMonitor;

    private AnalyzeGeneSets() {
    }

    /**
     * Command Line Tool to Analyze Gene Sets.
     *
     * @param args Command Line Arguments.  None Expected.
     * @throws java.io.IOException IO Error.
     */
    public static void main(String[] args) throws IOException {
        pMonitor = ProgressMonitor.getInstance();
        pMonitor.setConsoleMode(true);
        pMonitor.setCurrentMessage("Welcome to NetBox.  Initializing Database.  "
                + "Please wait a few moments...");

        System.out.println("Reading in Gene Sets");
        ReadGeneSets readGeneSets = new ReadGeneSets(new File("../data/c2.cp.v2.5.symbols.gmt"));
        //ReadGeneSets readGeneSets = new ReadGeneSets(new File("../data/OV_pathway_gene_sets.txt"));
        ArrayList<GeneSet> geneSetList = readGeneSets.getGeneSetList();
        System.out.println("Number of Gene Sets:  " + geneSetList.size());

        ArrayList<ByteProfileData> pList = new ArrayList<ByteProfileData>();
        CaseSetReader caseReader1 = new CaseSetReader(new File("../ova/cases_platsens.txt"));
        HashSet<String> platSensSet = caseReader1.getCaseIdSet();

        CaseSetReader caseReader2 = new CaseSetReader(new File("../ova/cases_platres.txt"));
        HashSet<String> platResistSet = caseReader2.getCaseIdSet();

        HashSet<String> unionCaseSet = new HashSet<String>();
        unionCaseSet.addAll(platSensSet);
        unionCaseSet.addAll(platResistSet);

        //  Get Mutation Data
        MutationReader mutReader = new MutationReader(new
                File("../ova/processed_final_mutations_6000.txt"),
                        unionCaseSet, true);
        ByteProfileData mutationData = mutReader.getByteProfile();
        pList.add(mutationData);

        //  Get Copy Number Data
        TabDelimReader cnaReader = new TabDelimReader(new File("../ova/data_CNA.txt"),
                GeneticAlterationType.COPY_NUMBER_ALTERATION, unionCaseSet, true);
        ByteProfileData cnaData = cnaReader.getByteProfile();
        pList.add(cnaData);

        //  Merge Profile Data
        ByteProfileMerger merger = new ByteProfileMerger(pList);
        ProfileDataSummary pSummary = new ProfileDataSummary(merger.getMergedProfile());

        //  Get Final Gene List
        ArrayList<GeneWithScore> geneWithScoreList = pSummary.getGeneFrequencyList();
        File alteredGenesFile = new File("altered_genes.txt");
        FileWriter writer = new FileWriter(alteredGenesFile);
        writer.write("FREQUENCY\n");
        ArrayList<String> geneList = new ArrayList<String>();
        for (GeneWithScore gene : geneWithScoreList) {
            geneList.add(gene.getGene());
            writer.write(gene.getGene() + "\t=" + gene.getScore() + "\n");
        }
        writer.close();


        ArrayList<String> platSensList = new ArrayList<String>();
        platSensList.addAll(platSensSet);

        ArrayList<String> platResistList = new ArrayList<String>();
        platResistList.addAll(platResistSet);

        writer = new FileWriter("gene_sets_platinum.txt");
        for (GeneSet geneSet : geneSetList) {
            ArrayList<Gene> genesInSet = geneSet.getGeneList();
            ArrayList<String> gList = new ArrayList<String>();
            StringBuffer geneBuf = new StringBuffer();
            for (Gene gene : genesInSet) {
                gList.add(gene.getGeneSymbol());
                geneBuf.append(gene.getGeneSymbol() + " ");
            }
            double sensitiveAffected = pSummary.getPercentCasesWhereGeneSetisAltered(gList,
                    platSensList);
            double resistantAffected = pSummary.getPercentCasesWhereGeneSetisAltered(gList,
                    platResistList);
            String str = geneSet.getName() + "\t" + geneBuf.toString() + "\t"
                    + sensitiveAffected + "\t" + resistantAffected
                    + "\n";
            writer.write(str);
            System.out.print(str);
        }
        writer.close();
        CommandLineUtil.showWarnings(pMonitor);
    }


}
