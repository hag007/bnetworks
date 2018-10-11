package org.mskcc.netbox.algorithm;

import edu.uci.ics.jung.graph.Graph;
import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.graph.InteractionToJung;
import org.mskcc.netbox.graph.InteractionUtil;
import org.mskcc.netbox.graph.NetworkPartitionState;
import org.mskcc.netbox.graph.NewmanGirvanModuleDetector;
import org.mskcc.netbox.model.Gene;
import org.mskcc.netbox.model.Interaction;
import org.mskcc.netbox.query.GeneQuery;
import org.mskcc.netbox.util.Formatter;
import org.mskcc.netbox.util.ProgressMonitor;
import org.mskcc.netbox.util.ReadGeneFile;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

/**
 * Random Null Model for NetBox.
 *
 * @author Ethan Cerami.
 */
public final class NetBoxRandom {
    private static final NumberFormat FORMATTER = Formatter.getDecimalFormat();
    private static final String TAB = "\t";
    private int modularityCounter = 0;
    private int sizeCounter = 0;

    private double observedNormalizedModularity;
    private int observedNetworkSize;
    private int numTrials;
    private ProgressMonitor pMonitor;

    /**
     * Creates N Random Networks and Evaluates their Modularity.
     *
     * @param file             Background Gene File.
     * @param numGenesPerTrial Number of Genes Selected Per Trial.
     * @param n                Total Number of Trials.
     * @param size             Observed Network Size.
     * @param q                Observed Network Modularity.
     * @throws IOException            IO Error.
     * @throws GraphCreationException Graph Creation Error.
     */
    public void createRandomNetworks(File file, Integer numGenesPerTrial, Integer n,
                                     int size, double q)
            throws IOException, GraphCreationException {
        Random randomGenerator = new Random();

        pMonitor = ProgressMonitor.getInstance();
        pMonitor.setCurrentMessage("Executing Background Model");
        HashMap<String, Gene> geneSymbolMap = GeneQuery.getGeneMapBySymbol();

        this.numTrials = n;
        this.observedNetworkSize = size;
        this.observedNormalizedModularity = q;

        FileWriter writer = new FileWriter(new File("random.dat"));

        //  Get background gene list.
        pMonitor.setCurrentMessage("Reading Background Gene File:  " + file.getAbsolutePath());
        ArrayList<String> backgroundGeneList = ReadGeneFile.readGeneFile(geneSymbolMap, file);
        pMonitor.setCurrentMessage("Background gene set has "
                + backgroundGeneList.size() + " genes.");
        for (int trial = 0; trial < numTrials; trial++) {
            pMonitor.setCurrentMessage("----------------------------");
            pMonitor.setCurrentMessage("Random Trial #" + trial);

            //  At each trial, randomly pick numGenesPerTrial genes from background
            ArrayList<String> currentGeneList = getRandomGeneSet(randomGenerator, numGenesPerTrial,
                    backgroundGeneList);

            StringBuffer geneBuf = new StringBuffer();
            for (String gene : currentGeneList) {
                geneBuf.append(gene + " ");
            }

            pMonitor.setCurrentMessage("The randomly selected genes are:  " + geneBuf.toString()
                    + " (Size:  " + currentGeneList.size() + ")");

            //  Connect these genes into a network
            ArrayList<Interaction> currentInteractionList = InteractionUtil.connectGenesNoLinkers(
                    currentGeneList);

            pMonitor.setCurrentMessage("Creates a network with " + currentInteractionList.size()
                    + " interactions.");

            //  Then partition the network into modules
            partitionNetwork(trial, writer, currentInteractionList);
        }
        writer.close();
        writeSummary();
    }

    private void writeSummary() throws IOException {
        pMonitor.setCurrentMessage("Summary of Results");
        pMonitor.setCurrentMessage("===================");

        File file = new File("summary.txt");
        FileWriter writer = new FileWriter(file);
        double pValue = sizeCounter / (double) numTrials;
        StringBuffer buf = new StringBuffer("Observed Network had size:  " + observedNetworkSize
                + " (p-value:  " + FORMATTER.format(pValue) + ").");
        writer.write(buf.toString() + "\n");
        pMonitor.setCurrentMessage(buf.toString());

        pValue = modularityCounter / (double) numTrials;
        buf = new StringBuffer("Observed Network had modularity:  "
                + FORMATTER.format(observedNormalizedModularity)
                + " (p-value:  " + FORMATTER.format(pValue) + ").");
        writer.write(buf.toString() + "\n");
        pMonitor.setCurrentMessage(buf.toString());

        buf = new StringBuffer("p-Values are based on " + numTrials + " random trials.");
        writer.write(buf.toString() + "\n");
        pMonitor.setCurrentMessage(buf.toString());

        writer.close();
        pMonitor.setCurrentMessage("Summary results written to:  " + file.getAbsolutePath());
    }

    private ArrayList<String> getRandomGeneSet(Random randomGenerator, int numGenesPerTrial,
                                               ArrayList<String> backgroundGeneList) {
        ArrayList<String> currentGeneList = new ArrayList<String>();
        for (int j = 0; j < numGenesPerTrial; j++) {
            int randomIndex = randomGenerator.nextInt(backgroundGeneList.size());
            currentGeneList.add(backgroundGeneList.get(randomIndex));
        }
        return currentGeneList;
    }

    private void partitionNetwork(int trialNum, FileWriter writer,
                                  ArrayList<Interaction> interactionList)
            throws GraphCreationException, IOException {
        if (interactionList.size() == 0) {
            writer.write(trialNum + TAB + "0" + TAB + "0\n");
            return;
        }
        Graph g = InteractionToJung.createGraph(interactionList);
        NewmanGirvanModuleDetector moduleDetector = new NewmanGirvanModuleDetector(g);
        NetworkPartitionState optimalState = moduleDetector.getOptimalPartitionState();
        pMonitor.setCurrentMessage("Max modularity occurs at:  "
                + optimalState.getNumEdgesRemoved()
                + " edge(s) removed.");
        pMonitor.setCurrentMessage("Results in:  " + optimalState.getNumModules()
                + " modules, with modularity of:  "
                + FORMATTER.format(optimalState.getNetworkModularity()));

        LocalRandomNullModel normalizedUtil = new LocalRandomNullModel(g,
                optimalState.getNetworkModularity(), 100);

        writer.write(trialNum + TAB + g.getVertices().size() + TAB
                + moduleDetector.isInitialStateRepresentedBySingleConnectedComponent()
                + TAB + FORMATTER.format(normalizedUtil.getNormalizedModularityScore())
                + "\n");

        double randomNetworkSize = g.getVertices().size();

        if (randomNetworkSize >= observedNetworkSize) {
            sizeCounter++;
        }
        if (normalizedUtil.getNormalizedModularityScore() > observedNormalizedModularity) {
            modularityCounter++;
        }

        writer.flush();
    }
}
