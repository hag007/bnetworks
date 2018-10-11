package org.mskcc.netbox.script;

import edu.uci.ics.jung.graph.Graph;
import org.mskcc.netbox.algorithm.GeneConnector;
import org.mskcc.netbox.algorithm.PercentLocalModuleDetector2;
import org.mskcc.netbox.genomic.ByteProfileData;
import org.mskcc.netbox.genomic.ByteProfileMerger;
import org.mskcc.netbox.genomic.GeneWithScore;
import org.mskcc.netbox.genomic.GeneticAlterationType;
import org.mskcc.netbox.genomic.ProfileDataSummary;
import org.mskcc.netbox.genomic.util.CaseSetReader;
import org.mskcc.netbox.genomic.util.MutationReader;
import org.mskcc.netbox.genomic.util.TabDelimReader;
import org.mskcc.netbox.report.HtmlReportGenerator;
import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.GlobalConfig;
import org.mskcc.netbox.util.GlobalSession;
import org.mskcc.netbox.util.NetworkWriter;
import org.mskcc.netbox.util.ProgressMonitor;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Main NetBox Command Line Program.
 *
 * @author Ethan Cerami.
 */
public final class NetBox2 {

    private NetBox2() {
    }

    /**
     * Command Line Main.
     *
     * @param args Command Line Arguments.
     * @throws Exception All Errors.
     */
    public static void main(String[] args) throws Exception {
        ProgressMonitor pMonitor = ProgressMonitor.getInstance();
        pMonitor.setConsoleMode(true);

        if (args.length == 0) {
            System.out.println("usage:  netBox2.py netbox.props");
            System.exit(1);
        }

        try {
            GlobalConfig config = GlobalConfig.getInstance();
            config.loadProperties(new File(args[0]));

            //  Make this call here to init the database
            pMonitor.setCurrentMessage("Welcome to NetBox.  Initializing Database.");
            GlobalSession globalSession = GlobalSession.getInstance();

            CaseSetReader caseReader = new CaseSetReader(config.getCaseFile());
            HashSet<String> caseIdSet = caseReader.getCaseIdSet();

            ArrayList<ByteProfileData> pList = new ArrayList<ByteProfileData>();

            MutationReader mutReader = new MutationReader(config.getMutationFile(),
                    caseIdSet, true);
            ByteProfileData mutationData = mutReader.getByteProfile();
            pList.add(mutationData);

            TabDelimReader cnaReader = new TabDelimReader(config.getCnaFile(),
                    GeneticAlterationType.COPY_NUMBER_ALTERATION, caseIdSet, true);
            ByteProfileData cnaData = cnaReader.getByteProfile();
            pList.add(cnaData);

            if (config.getMRNAFile() != null) {
                TabDelimReader mrnaReader = new TabDelimReader(config.getMRNAFile(),
                        GeneticAlterationType.MRNA_EXPRESSION, caseIdSet, true);
                pList.add(mrnaReader.getByteProfile());
            }

            CommandLineUtil.showWarnings(pMonitor);
            ByteProfileMerger merger = new ByteProfileMerger(pList);

            ProfileDataSummary pSummary = new ProfileDataSummary(merger.getMergedProfile());

            ArrayList<GeneWithScore> geneWithScoreList = pSummary.getGeneFrequencyList();
            FileWriter writer = new FileWriter("altered_genes.txt");
            ArrayList<String> geneList = new ArrayList<String>();
            for (GeneWithScore gene : geneWithScoreList) {
                if (gene.getScore() > config.getGeneFrequencyThreshold()) {
                    geneList.add(gene.getGene());
                    writer.write(gene.getGene() + "\t" + gene.getScore() + "\n");
                }
            }
            writer.close();

            GeneConnector geneConnector = new GeneConnector(geneList,
                    config.getShortestPathThreshold(), config.getPValueCutOff());
            Graph g = geneConnector.getGraph();

            PercentLocalModuleDetector2 moduleDetector = new PercentLocalModuleDetector2(g,
                    pSummary, config.getSeedRadiusThreshold(), config.getDeltaThreshold(),
                    config.getModuleFrequencyThreshold(), config.getNumLocalTrials());

            pMonitor.setCurrentMessage("------------------------------------");
            NetworkWriter.outputNetwork(g);
            NetworkWriter.outputModulesByPercentage(moduleDetector);
            NetworkWriter.outputGeneAttributes(geneWithScoreList);
            NetworkWriter.outputNodeAttributes(geneConnector.getAlteredGeneList(),
                    geneConnector.getLinkerList());

            HtmlReportGenerator htmlGenerator = HtmlReportGenerator.getInstance();
            htmlGenerator.appendGraph(geneList, g, geneConnector);
            htmlGenerator.appendPercentModuleDetection(moduleDetector);
            htmlGenerator.appendLinkers(geneConnector.getLinkerList());
            htmlGenerator.finalizeReport();
        } catch (Throwable t) {
            exitCommandLine(t);
        }
    }

    /**
     * Exits the Command Line with the Specified Error Message.
     *
     * @param e Exception Object.
     */
    private static void exitCommandLine(Throwable e) {
        System.out.println("Error:  " + e.getMessage());
        e.printStackTrace();
        System.exit(1);
    }

}
