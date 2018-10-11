package org.mskcc.netbox.algorithm;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import org.apache.log4j.Logger;
import org.mskcc.netbox.genomic.GeneWithScore;
import org.mskcc.netbox.genomic.ProfileDataSummary;
import org.mskcc.netbox.graph.Module;
import org.mskcc.netbox.stats.BenjaminiHochbergFDR;
import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.ProgressMonitor;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;


/**
 * Detects Modules, based on Percent Alteration Values.
 *
 * @author Ethan Cerami.
 */
public final class PercentLocalModuleDetector2 implements ModuleDetector {
    private Graph graph;
    private ProfileDataSummary pSummary;
    private StringLabeller labeller;
    private int maxD;
    private ArrayList<Module> moduleList = new ArrayList<Module>();
    private ProgressMonitor pMonitor = ProgressMonitor.getInstance();
    private static Logger logger = Logger.getLogger(PercentLocalModuleDetector.class);
    private double moduleFrequencyThreshold;
    private double deltaThreshold;
    private int numRandomTrials;

    /**
     * Constructor.
     *
     * @param g                   Graph of Interest.
     * @param p                   Profile Data Summary Object.
     * @param d                   Max Distance from Seed.
     * @param dThreshold          Delta Threshold for Adding New Nodes to SubNet.
     * @param mFrequencyThreshold Module Frequency Threshold.
     * @param n                   Number of Random Trials.
     */
    public PercentLocalModuleDetector2(Graph g, ProfileDataSummary p, int d, double dThreshold,
                                       double mFrequencyThreshold, int n) {
        graph = g;
        labeller = StringLabeller.getLabeller(graph);
        pSummary = p;
        maxD = d;
        deltaThreshold = dThreshold;
        moduleFrequencyThreshold = mFrequencyThreshold;
        numRandomTrials = n;
        execute();
    }

    /**
     * Gets Discovered Modules.
     *
     * @return ArrayList of Modules.
     */
    public ArrayList<Module> getModules() {
        return moduleList;
    }

    private void execute() {

        int currentModuleId = 0;
        ArrayList<GeneWithScore> geneScoreList = pSummary.getGeneFrequencyList();
        HashMap<String, GeneWithScore> geneScoreMap = new HashMap<String, GeneWithScore>();
        for (GeneWithScore geneScore : geneScoreList) {
            geneScoreMap.put(geneScore.getGene(), geneScore);
        }


        //  Iterate through all genes, beginning with most frequenly altered.

        pMonitor.setCurrentMessage("Examining All Subnetworks");
        pMonitor.setMaxValue(graph.numVertices());
        Set<Vertex> vertexSet = graph.getVertices();

        //  Examine each local neighborhood
        for (Vertex vertex : vertexSet) {

            String label = labeller.getLabel(vertex);

            System.out.println("Examining neighborhood centered at:  " + label);

            ArrayList<GeneWithScore> genesInModuleList = new ArrayList<GeneWithScore>();

            GeneWithScore currentGene = new GeneWithScore();
            currentGene.setGene(label);

            ArrayList<String> subNet = new ArrayList<String>();
            subNet.add(currentGene.getGene());

            Set<Vertex> neighborSet = vertex.getNeighbors();
            for (Vertex v : neighborSet) {
                String s = labeller.getLabel(v);
                subNet.add(s);
            }


            double frequency = pSummary.getPercentCasesWhereGeneSetisAltered(subNet);

            HashSet<String> genesExploredSet = new HashSet<String>();
            for (String gene : subNet) {
                GeneWithScore tempGene = new GeneWithScore();
                tempGene.setGene(gene);
                tempGene.setScore(frequency);
                genesInModuleList.add(tempGene);
                genesExploredSet.add(gene);
            }

            //  Only Store Modules of Size > 1
            if (frequency > moduleFrequencyThreshold) {
                Module module = new Module();
                module.setModuleId(currentModuleId);
                module.setGeneList(genesInModuleList);
                module.setScore(frequency);
                module.setGenesExplored(genesExploredSet);
                moduleList.add(module);
                currentModuleId++;
            }
            pMonitor.incrementCurValue();
            CommandLineUtil.showProgress(pMonitor);
        }
        executeRandomBackgroundModel();
    }

    private void executeRandomBackgroundModel() {
        Random randomGenerator = new Random();
        ArrayList<String> observedGeneList = pSummary.getObservedGeneList();

        pMonitor.setCurrentMessage("Executing Random Background Model");
        pMonitor.setMaxValue(moduleList.size());

        for (Module module : moduleList) {

            int counter = 0;
            for (int i = 0; i < numRandomTrials; i++) {

                //  Get Observed Frequency Score
                double observedScore = module.getScore();

                //  Pick random gene set of size N.
                int geneSetSize = module.getGenesExplored().size() + 1;
                ArrayList<String> randomGeneList = new ArrayList<String>();
                for (int j = 0; j < geneSetSize; j++) {
                    int index = randomGenerator.nextInt(observedGeneList.size());
                    randomGeneList.add(observedGeneList.get(index));
                }

                double randomScore = pSummary.getPercentCasesWhereGeneSetisAltered(randomGeneList);
                if (randomScore >= observedScore) {
                    counter++;
                }
            }
            double pValue = counter / (double) numRandomTrials;
            module.setPValueUnAdjusted(pValue);
            logger.info("Module:  " + module.getModuleId() + ", "
                    + module.getLabel() + ", frequency:  " + module.getScore() + ", "
                    + " p-value:  " + module.getPValueUnAdjusted());
            pMonitor.incrementCurValue();
            CommandLineUtil.showProgress(pMonitor);
        }

        //  Sort all Modules by P-Value
        Collections.sort(moduleList, new ModuleComparator());

        double[] pValues = new double[moduleList.size()];
        for (int i = 0; i < moduleList.size(); i++) {
            Module module = moduleList.get(i);
            pValues[i] = module.getPValueUnAdjusted();
        }

        //  Adjust for multiple hypothesis via Benjamani Hochberg FDR
        pMonitor.setCurrentMessage("Calculating Benjamini Hochberg FDR Adjustment");
        BenjaminiHochbergFDR fdr = new BenjaminiHochbergFDR(pValues);
        fdr.calculate();
        double[] adjustedPValues = fdr.getAdjustedPvalues();

        //  Reset the IDs.
        for (int i = 0; i < moduleList.size(); i++) {
            Module module = moduleList.get(i);
            module.setModuleId(i);
            module.setPValueFdrAdjusted(adjustedPValues[i]);
        }
    }

    /**
     * Returns the Gene that Results in the Best Frequency of Alteration Improvement,
     * or null if none of them do.
     *
     * @param subNet      Current SubNet.
     * @param neighborSet Current Neighbor Set.
     * @return Best Addition.
     */
    private GeneWithScore getBestAddition(HashSet<String> subNet, HashSet<String> neighborSet) {

        //  First, determine base line frequency score.
        ArrayList<String> subListList = convertHashSetToArrayList(subNet);
        double baseScore = pSummary.getPercentCasesWhereGeneSetisAltered(subListList);
        logger.info("Base line frequency is:  " + baseScore);

        //  Evaluate all neighbor candidates
        double maxScore = baseScore;
        String maxCandidate = null;
        for (String candidate : neighborSet) {
            ArrayList<String> neighborList = convertHashSetToArrayList(subNet);
            neighborList.add(candidate);
            double candidateFrequency = pSummary.getPercentCasesWhereGeneSetisAltered(neighborList);
            logger.info("Adding gene:  " + candidate + " would result in frequency:  "
                    + candidateFrequency);
            if (candidateFrequency > maxScore) {
                maxScore = candidateFrequency;
                maxCandidate = candidate;
            }
        }
        if (maxCandidate != null && (maxScore - baseScore) > deltaThreshold) {
            GeneWithScore geneWithScore = new GeneWithScore();
            geneWithScore.setGene(maxCandidate);
            geneWithScore.setScore(maxScore);
            logger.info("Best candidate to add is:  " + geneWithScore);
            return geneWithScore;
        } else {
            logger.info("No good candidates to add.  Module exploration at end.");
            return null;
        }
    }

    /**
     * Converts a HashSet to an ArrayList.
     *
     * @param set HashSet.
     * @return ArrayList.
     */
    private ArrayList<String> convertHashSetToArrayList(HashSet<String> set) {
        ArrayList<String> list = new ArrayList<String>();
        for (String s : set) {
            list.add(s);
        }
        return list;
    }
}
