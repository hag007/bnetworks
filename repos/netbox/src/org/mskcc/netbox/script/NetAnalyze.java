package org.mskcc.netbox.script;

import edu.uci.ics.jung.graph.Graph;
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

/**
 * Main GenesToNetwork Command Line Program.
 *
 * @author Ethan Cerami.
 */
public final class NetAnalyze {
    private static NumberFormat formatter = Formatter.getDecimalFormat();
    private static String cmdLineUsage = "command line usage:  netAnalyze.py netbox.props";
    private ArrayList<String> geneList;
    private ProgressMonitor pMonitor;
    private VertexEdgePair largestComponentInfo;

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

            NetAnalyze netAnalyze = new NetAnalyze();
            netAnalyze.execute();
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
        pMonitor.setCurrentMessage("Welcome to NetBox.  Initializing Database.  "
                + "Please wait a few moments...");
        GlobalSession globalSession = GlobalSession.getInstance();

        //  Get Global Config Options
        GlobalConfig globalConfig = GlobalConfig.getInstance();
        File geneFile = globalConfig.getGeneFile();
        int shortestPathThreshold = globalConfig.getShortestPathThreshold();
        double pValueCutOff = globalConfig.getPValueCutOff();

        //  Get Maps of Genes
        HashMap<String, Gene> geneSymbolMap = GeneQuery.getGeneMapBySymbol();

        File alteredGenesFile = null;

        if (geneFile != null) {
            //  Connect Genes
            geneList = ReadGeneFile.readGeneFile(geneSymbolMap, geneFile);
            pMonitor.setCurrentMessage("Total number of genes read in:  "
                    + geneList.size() + ".");
        } else {
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
        }

        GeneConnector geneConnector;
        geneConnector = new GeneConnector(geneList, shortestPathThreshold,
                pValueCutOff);
        if (globalConfig.isInteractiveMode()) {
            System.out.print("\nEnter [1] to accept network;  [2] to enter new p-value "
                    + "threshold:  ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            String line = br.readLine();
            while (line.equals("2")) {
                System.out.print("Enter new p-value threshold:  ");
                line = br.readLine();
                double p = Double.parseDouble(line);
                globalConfig.setPValueCutOff(p);
                geneConnector = new GeneConnector(geneList, 2, p);
                System.out.print("\nEnter [1] to accept network;  [2] to enter new p-value "
                        + "threshold:  ");
                line = br.readLine();
            }
        }

        //  Partition Graph into Modules via chosen algorithm
        Graph g = geneConnector.getGraph();
        HtmlReportGenerator html = HtmlReportGenerator.getInstance();
        html.appendGraph(geneList, g, geneConnector);

        String algo = globalConfig.getNetworkPartitionAlgorithm();
        double observedNetworkModularity = 0;
        if (globalConfig.identifyModules()) {
            if (algo.equals(GlobalConfig.NG)) {
                observedNetworkModularity = executeNewmanGirvan(geneConnector);
            } else {
                observedNetworkModularity = executeSA(geneConnector);
            }
        }

        if (g.getVertices().size() > 0) {
            NetworkWriter.outputNetwork(g);
            NetworkWriter.outputNodeAttributes(geneConnector.getAlteredGeneList(),
                    geneConnector.getLinkerList());
        }

        if (globalConfig.identifyModules()) {
            largestComponentInfo = LargestComponentUtil.determineSizeOfLargestComponent(g);
        }

        if (alteredGenesFile != null) {
            pMonitor.setCurrentMessage("Attribute file containing gene alteration frequencies, "
                    + "suitable for loading into Cytoscape is available at:  "
                    + alteredGenesFile.getAbsolutePath());
        }

        //  Execute Local Null Model
        if (globalConfig.identifyModules()) {
            int numLocalTrials = globalConfig.getNumLocalTrials();
            if (numLocalTrials > 0) {
                StringBuffer buf = new StringBuffer();
                pMonitor.setCurrentMessage("Executing Local Random Null Model");
                LocalRandomNullModel nullModel = new LocalRandomNullModel(g, observedNetworkModularity,
                        numLocalTrials);
                buf.append("Observed Network Modularity is:  " + observedNetworkModularity + "\n");
                buf.append("Scaled Network Modularity is:  "
                        + nullModel.getNormalizedModularityScore() + "\n");
                buf.append("Random Networks:  mean modularity = " + nullModel.getRandomMean()
                        + ", sd = " + nullModel.getRandomStandardDeviation() + "\n");
                buf.append("Based on " + numLocalTrials + " random trials." + "\n");
                buf.append("Network modularity for all locally rewired networks available at:  "
                        + "<a href='" + nullModel.getFile().getAbsolutePath() + "'>"
                        + nullModel.getFile().getAbsolutePath() + "</a>");
                html.appendLocalNullModel(buf.toString());
            }

            //  Execute Global Null Model
            int numGlobalTrials = globalConfig.getNumGlobalTrials();
            if (numGlobalTrials > 0) {
                StringBuffer buf = new StringBuffer();
                pMonitor.setCurrentMessage("Executing Global Random Null Model");
                GlobalRandomNullModel nullModel = new GlobalRandomNullModel(geneConnector,
                        numGlobalTrials, shortestPathThreshold, pValueCutOff, largestComponentInfo);
                buf.append("Observed network has largest component with "
                        + largestComponentInfo.getNumVertices() + " vertices "
                        + "(p-value:  " + nullModel.getPValueNodes() + ")");
                buf.append(" and " + largestComponentInfo.getNumEdges()
                        + " edges (p-value:  " + nullModel.getPValueEdges() + ").\n");
                buf.append("Based on " + numGlobalTrials + " random trials.\n");
                buf.append("Number of vertices / edges in the largest component in all random "
                        + "networks available at:  <a href='" + nullModel.getFile().getAbsolutePath()
                        + "'>" + nullModel.getFile().getAbsolutePath() + "</a>");
                html.appendGlobalNullModel(buf.toString());
            }
        }

        html.appendLinkers(geneConnector.getLinkerList());
        html.finalizeReport();
    }


    private double executeSA(GeneConnector geneConnector) throws GraphCreationException,
            IOException {
        HtmlReportGenerator html = HtmlReportGenerator.getInstance();
        NetCartoAnnealing sa = new NetCartoAnnealing(geneConnector.getGraph());
        sa.execute();
        NetworkWriter.outputModules(sa);
        html.appendModuleDetection(sa, geneConnector.getLinkerList());
        return sa.getFinalModularity();
    }

    /**
     * Parition the Graph into Modules via Newman Girvan Algorithm.
     *
     * @param geneConnector GeneConnector Object.
     * @throws IOException IO Error.
     */
    private double executeNewmanGirvan(GeneConnector geneConnector) throws IOException,
            GraphCreationException {
        Graph g = geneConnector.getGraph();

        pMonitor.setCurrentMessage("\nBegin Module Detection (" + "Number of nodes:  "
                + g.getVertices().size()
                + ", Number of edges:  " + g.getEdges().size() + ")");

        //  Create Module Detector
        NewmanGirvanModuleDetector moduleDetector = new NewmanGirvanModuleDetector(g);

        //  Get Partitioning that Results in Optimal Modularity
        NetworkPartitionState optimalState = moduleDetector.getOptimalPartitionState();
        pMonitor.setCurrentMessage("Max modularity occurs at:  "
                + optimalState.getNumEdgesRemoved()
                + " edge(s) removed.");
        pMonitor.setCurrentMessage("Results in:  " + optimalState.getNumModules()
                + " modules, with modularity of:  "
                + formatter.format(optimalState.getNetworkModularity()));

        NetworkWriter.outputModules(moduleDetector);
        HtmlReportGenerator html = HtmlReportGenerator.getInstance();
        html.appendModuleDetection(moduleDetector, geneConnector.getLinkerList());
        return optimalState.getNetworkModularity();
    }
}

/**
 * NetBox Command Line Options.
 */
class NetBoxOptions extends OptionParser {

    /**
     * Constructor.
     */
    NetBoxOptions() {
        OptionSpecBuilder builder = accepts("d", "Prints diagnostics in the event of an error.  "
                + "Will print error stack trace.");
        builder.withOptionalArg();

        accepts("i", "Interactive mode.");

        accepts("h", "Show help.");
    }
}

