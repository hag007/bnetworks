package org.mskcc.netbox.script;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import joptsimple.OptionException;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpecBuilder;
import org.mskcc.netbox.algorithm.GeneConnector;
import org.mskcc.netbox.algorithm.GlobalRandomNullModel;
import org.mskcc.netbox.algorithm.LargestComponentUtil;
import org.mskcc.netbox.algorithm.LocalRandomNullModel;
import org.mskcc.netbox.algorithm.VertexEdgePair;
import org.mskcc.netbox.genomic.ByteProfileData;
import org.mskcc.netbox.genomic.ByteProfileMerger;
import org.mskcc.netbox.genomic.GeneWithScore;
import org.mskcc.netbox.genomic.GeneticAlterationType;
import org.mskcc.netbox.genomic.ProfileDataSummary;
import org.mskcc.netbox.genomic.util.CaseSetReader;
import org.mskcc.netbox.genomic.util.MutationReader;
import org.mskcc.netbox.genomic.util.TabDelimReader;
import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.graph.NetworkPartitionState;
import org.mskcc.netbox.graph.NewmanGirvanModuleDetector;
import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.netcarto.NetCartoAnnealing;
import org.mskcc.netbox.query.GeneQuery;
import org.mskcc.netbox.report.HtmlReportGenerator;
import org.mskcc.netbox.util.Formatter;
import org.mskcc.netbox.util.GlobalConfig;
import org.mskcc.netbox.util.GlobalSession;
import org.mskcc.netbox.util.NetworkWriter;
import org.mskcc.netbox.util.ProgressMonitor;
import org.mskcc.netbox.util.ReadGeneFile;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import jsc.correlation.KendallCorrelation;
import jsc.datastructures.PairedData;

/**
 * Main BackBone Command Line Program.
 *
 * @author Ethan Cerami.
 */
public final class BackBone {
    private static String cmdLineUsage = "command line usage:  backBone.py netbox.props";
    private ArrayList<String> geneList;
    private ProgressMonitor pMonitor;

    /**
     * Command Line Main.
     *
     * @param args Command Line Arguments.
     */
    public static void main(String[] args) {
        OptionParser parser = new NetBoxOptions();
        GlobalConfig globalConfig = GlobalConfig.getInstance();
        GlobalSession globalSession = GlobalSession.getInstance();

        OptionSet options = null;
        try {
            if (args.length == 0) {
                System.out.println("\n" + cmdLineUsage + "\n");
                parser.printHelpOn(System.out);
                System.exit(1);
            }

            options = parser.parse(args);
            if (options.has("h")) {
                parser.printHelpOn(System.out);
                System.exit(1);
            }
            if (options.has("d")) {
                globalConfig.setDebugMode(true);
            }
            if (options.has("i")) {
                globalConfig.setInteractiveMode(true);
            }
            globalConfig.loadProperties(new File(args[0]));

            BackBone backBone = new BackBone();
            backBone.execute();
        } catch (OptionException e) {
            exitCommandLine(e);
        } catch (IllegalArgumentException e) {
            exitCommandLine(e);
        } catch (IOException e) {
            exitCommandLine(e);
        } catch (Throwable t) {
            exitCommandLine(t);
        } finally {
            globalSession.closeAll();
        }

    }

    /**
     * Exits the Command Line with the Specified Error Message.
     *
     * @param e Exception Object.
     */
    private static void exitCommandLine(Throwable e) {
        System.out.println("Command Line Error:  " + e.getMessage() + "  Use -h to get help.");
        if (GlobalConfig.getInstance().isDebugMode()) {
            e.printStackTrace();
        }
        System.exit(1);
    }

    /**
     * Executes NetBox.
     *
     * @throws GraphCreationException Graph Creation Error.
     * @throws IOException            IO Error.
     */
    public void execute() throws GraphCreationException, IOException {
        pMonitor = ProgressMonitor.getInstance();
        pMonitor.setConsoleMode(true);
        pMonitor.setCurrentMessage("Welcome to NetBox BackBone.  Initializing Database.  "
                + "Please wait a few moments...");
        GlobalSession globalSession = GlobalSession.getInstance();

        //  Get Global Config Options
        GlobalConfig globalConfig = GlobalConfig.getInstance();
        int shortestPathThreshold = globalConfig.getShortestPathThreshold();
        double pValueCutOff = globalConfig.getPValueCutOff();

        File alteredGenesFile = null;

        ArrayList<ByteProfileData> pList = new ArrayList<ByteProfileData>();
        CaseSetReader caseReader = new CaseSetReader(globalConfig.getCaseFile());
        HashSet<String> caseIdSet = caseReader.getCaseIdSet();

        //  Get Mutation Data
        MutationReader mutReader = new MutationReader(globalConfig.getMutationFile(),
                caseIdSet, true);
        ByteProfileData mutationData = mutReader.getByteProfile();
        pList.add(mutationData);

        //  Get Copy Number Data
        if (globalConfig.getCnaFile() != null) {
            TabDelimReader cnaReader = new TabDelimReader(globalConfig.getCnaFile(),
                    GeneticAlterationType.COPY_NUMBER_ALTERATION, caseIdSet, true);
            ByteProfileData cnaData = cnaReader.getByteProfile();
            pList.add(cnaData);
        }

        //  Merge Profile Data
        ByteProfileMerger merger = new ByteProfileMerger(pList);
        ProfileDataSummary pSummary = new ProfileDataSummary(merger.getMergedProfile());

        //  Get Final Gene List
        ArrayList<GeneWithScore> geneWithScoreList = pSummary.getGeneFrequencyList();
        alteredGenesFile = new File("altered_genes.txt");
        FileWriter writer = new FileWriter(alteredGenesFile);
        writer.write("FREQUENCY\n");
        geneList = new ArrayList<String>();
        for (GeneWithScore gene : geneWithScoreList) {
            if (gene.getScore() > globalConfig.getGeneFrequencyThreshold()) {
                geneList.add(gene.getGene());
                writer.write(gene.getGene() + "\t=" + gene.getScore() + "\n");
            }
        }
        writer.close();

        GeneConnector geneConnector;
        geneConnector = new GeneConnector(geneList, shortestPathThreshold,
                pValueCutOff);

        Graph g = geneConnector.getGraph();
        traceBackBone(g, pSummary);

        HtmlReportGenerator html = HtmlReportGenerator.getInstance();
        html.appendGraph(geneList, g, geneConnector);

        if (alteredGenesFile != null) {
            pMonitor.setCurrentMessage("Attribute file containing gene alteration frequencies, "
                    + "suitable for loading into Cytoscape is available at:  "
                    + alteredGenesFile.getAbsolutePath());
        }

        html.appendLinkers(geneConnector.getLinkerList());
        html.finalizeReport();
    }

    private void traceBackBone (Graph g, ProfileDataSummary pSummary) throws IOException {
        HashSet<String> edgeSet = new HashSet<String>();
        FileWriter writer1 = new FileWriter ("network.sif");
        FileWriter writer2 = new FileWriter ("edges.txt");
        writer2.write("KENDALL\n");
        StringLabeller labeller = StringLabeller.getLabeller(g);
        Set<Vertex> vertexSet = g.getVertices();
        for (Vertex vertex:  vertexSet) {
            String geneA = labeller.getLabel(vertex);
            Set<Vertex> neighborSet = vertex.getNeighbors();
            for (Vertex neighbor:  neighborSet) {
                int n=0;
                String geneB = labeller.getLabel(neighbor);
                ArrayList<String> caseList = pSummary.getObservedCaseList();
                ArrayList<Double> xList = new ArrayList<Double>();
                ArrayList<Double> yList = new ArrayList<Double>();
                for (String caseId:  caseList) {
                    boolean geneAAltered = pSummary.isGeneAltered(geneA, caseId);
                    boolean geneBAltered = pSummary.isGeneAltered(geneB, caseId);
                    if (geneAAltered || geneBAltered) {
                        if (geneAAltered) {
                            xList.add(1.0);
                        } else {
                            xList.add(0.0);
                        }
                        if (geneBAltered) {
                            yList.add(1.0);
                        } else {
                            yList.add(0.0);
                        }
                        n++;
                    }
                }

                double x[] = new double[n];
                double y[] = new double[n];
                for (int i=0; i<xList.size(); i++) {
                    x[i] = xList.get(i);
                    y[i] = yList.get(i);
                }

                PairedData pairedData = new PairedData(x,y);
                KendallCorrelation kendall = new KendallCorrelation(pairedData);
                if (kendall.getR() < -.7) {
                    String key = getKey(geneA, geneB);
                    if (!edgeSet.contains(key)) {
                        writer1.write (geneA + " pp " + geneB + "\n");
                        writer2.write (geneA + " (pp) " + geneB + "= " + kendall.getR() + "\n");
                        edgeSet.add(key);
                    }
                }
            }
        }
        writer1.close();
        writer2.close();
    }

    private String getKey (String geneA, String geneB) {
        if (geneA.compareTo(geneB) > 0) {
            return geneA + ":" + geneB;
        } else {
            return geneB + ":" + geneA;
        }
    }
}

/**
 * BackBone Command Line Options.
 */
class BackBoneOptions extends OptionParser {

    /**
     * Constructor.
     */
    BackBoneOptions() {
        OptionSpecBuilder builder = accepts("d", "Prints diagnostics in the event of an error.  "
                + "Will print error stack trace.");
        builder.withOptionalArg();

        accepts("i", "Interactive mode.");

        accepts("h", "Show help.");
    }
}

