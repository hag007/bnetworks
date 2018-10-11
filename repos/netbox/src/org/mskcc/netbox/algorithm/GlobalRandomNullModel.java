package org.mskcc.netbox.algorithm;

import org.mskcc.netbox.graph.GraphCreationException;
import org.mskcc.netbox.graph.NetworkStatsUtil;
import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.GlobalConfig;
import org.mskcc.netbox.util.ProgressMonitor;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

/**
 * Global Random Null Model.
 *
 * @author Ethan Cerami.
 */
public final class GlobalRandomNullModel {
    private double pValueNodes;
    private double pValueEdges;
    private File file;

    /**
     * Constructor.
     *
     * @param geneConnector     Gene Connector Object.
     * @param numTrials         Number of Trials to Execute.
     * @param shortestPathThreshold Shortest Path Threshold.
     * @param pValueCutOff      P-Value Threshold.
     * @param largestComponentInfo  Observed Largest Component Size.
     * @throws GraphCreationException   Graph Creation Error.
     * @throws IOException  IO Error.
     */
    public GlobalRandomNullModel(GeneConnector geneConnector, int numTrials,
                                 int shortestPathThreshold, double pValueCutOff,
                                 VertexEdgePair largestComponentInfo)
            throws GraphCreationException, IOException {
        ProgressMonitor pMonitor = ProgressMonitor.getInstance();
        Random randomGenerator = new Random();
        int size = geneConnector.getAlteredGeneList().size();
        ArrayList<String> genesInNetwork = NetworkStatsUtil.getInstance().getGenesInNetwork();

        file = new File("global_null.txt");
        FileWriter writer = new FileWriter(file);

        GlobalConfig config = GlobalConfig.getInstance();
        boolean originalQuiet = config.isVeryQuiet();
        config.setBeVeryQuiet(true);

        int nodeCounter = 0;
        int edgeCounter = 0;
        pMonitor.setCurrentMessage("Executing Global Null Model");
        pMonitor.setMaxValue(numTrials);
        for (int i = 0; i < numTrials; i++) {
            ArrayList<String> randomGeneList = new ArrayList<String>();
            for (int j = 0; j < size; j++) {
                int randomIndex = randomGenerator.nextInt(genesInNetwork.size());
                String gene = genesInNetwork.get(randomIndex);
                randomGeneList.add(gene);
            }
            GeneConnector geneConnector2 = new GeneConnector(randomGeneList,
                    shortestPathThreshold, pValueCutOff);
            VertexEdgePair randomPair =
                    LargestComponentUtil.determineSizeOfLargestComponent(
                            geneConnector2.getGraph());
            if (randomPair.getNumVertices() >= largestComponentInfo.getNumVertices()) {
                nodeCounter++;
            }
            if (randomPair.getNumEdges() >= largestComponentInfo.getNumEdges()) {
                edgeCounter++;
            }
            writer.write(randomPair.getNumVertices() + "\t" + randomPair.getNumEdges() + "\n");
            writer.flush();
            pMonitor.incrementCurValue();
            pMonitor.setCurrentMessage("Executing Global Null Model");
            CommandLineUtil.showProgress(pMonitor);
        }
        pValueNodes = nodeCounter / (double) numTrials;
        pValueEdges = edgeCounter / (double) numTrials;
        writer.close();

        config.setBeVeryQuiet(originalQuiet);
    }

    /**
     * Gets File containing null model results.
     * @return file.
     */
    public File getFile() {
        return file;
    }

    /**
     * Gets P-Value Associated with Node Size.
     * @return p-value.
     */
    public double getPValueNodes() {
        return pValueNodes;
    }

    /**
     * Gets P-Value Associated with Edge Size.
     * @return p-value.
     */
    public double getPValueEdges() {
        return pValueEdges;
    }
}
