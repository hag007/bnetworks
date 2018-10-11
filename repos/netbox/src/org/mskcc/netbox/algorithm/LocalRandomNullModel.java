package org.mskcc.netbox.algorithm;

import edu.uci.ics.jung.graph.Graph;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.mskcc.netbox.graph.GraphPermute;
import org.mskcc.netbox.graph.NetworkPartitionState;
import org.mskcc.netbox.graph.NewmanGirvanModuleDetector;
import org.mskcc.netbox.util.Formatter;
import org.mskcc.netbox.util.GlobalConfig;
import org.mskcc.netbox.util.ProgressMonitor;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;

/**
 * Utility Class for Creating a Normalized Z-Score for Network Modularity.
 *
 * @author Ethan Cerami.
 */
public final class LocalRandomNullModel {
    private static NumberFormat formatter = Formatter.getDecimalFormat();
    private double zScore;
    private StandardDeviation sd;
    private double randomMean = 0;
    private File file;

    /**
     * Constructor.
     *
     * @param g         Graph Original Graph.
     * @param qObserved Observed Modularity Score.
     * @param numTrials Number of Random Trials.
     * @throws IOException IO Error.
     */
    public LocalRandomNullModel(Graph g, double qObserved, int numTrials)
            throws IOException {

        sd = new StandardDeviation();
        file = new File("local_null.txt");
        FileWriter writer = new FileWriter(file);
        GraphPermute graphPermute = new GraphPermute(g);

        // Perform N Random Trials
        for (int i = 0; i < numTrials; i++) {

            //  At Each Random Trial, Permute the Local Wiring of the Graph
            Graph randomGraph = graphPermute.next();
            ProgressMonitor.getInstance().setCurrentMessage("Random Network #" + i
                    + " (Number of Nodes:  " + randomGraph.getVertices().size()
                    + ", Number of Edges:  " + randomGraph.getEdges().size() + ").");

            //  And, then calculate network modularity
            GlobalConfig config = GlobalConfig.getInstance();
            boolean originalQuiet = config.isVeryQuiet();
            config.setBeVeryQuiet(true);

            NewmanGirvanModuleDetector moduleDetector = new NewmanGirvanModuleDetector(randomGraph);
            NetworkPartitionState optimalState = moduleDetector.getOptimalPartitionState();

            config.setBeVeryQuiet(originalQuiet);
            ProgressMonitor.getInstance().setCurrentMessage("Random Network #"
                    + i + ": unscaled modularity score:  "
                    + formatter.format(optimalState.getNetworkModularity()));
            writer.write(optimalState.getNetworkModularity() + "\n");
            writer.flush();

            randomMean += optimalState.getNetworkModularity();
            sd.increment(optimalState.getNetworkModularity());
        }
        writer.close();

        //  Finally, calculate standard Z-Score
        randomMean = randomMean / numTrials;
        zScore = (qObserved - randomMean) / sd.getResult();
    }

    /**
     * Gets the File containing Results.
     * @return file.
     */
    public File getFile() {
        return file;
    }

    /**
     * Gets the Normalized Modularity Score.
     *
     * @return normalized modularity score.
     */
    public double getNormalizedModularityScore() {
        return zScore;
    }

    /**
     * Gets the Random Mean.
     *
     * @return Random Mean.
     */
    public double getRandomMean() {
        return randomMean;
    }

    /**
     * Gets the Random Standard Deviation.
     *
     * @return Random Standard Deviation.
     */
    public double getRandomStandardDeviation() {
        return sd.getResult();
    }
}
