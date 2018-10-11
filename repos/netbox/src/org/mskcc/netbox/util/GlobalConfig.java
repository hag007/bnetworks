package org.mskcc.netbox.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

/**
 * Global Configuration Object.
 *
 * @author Ethan Cerami
 */
public final class GlobalConfig {
    /**
     * Static Global Config Singleton.
     */
    private static GlobalConfig globalConfig;

    private static final int DEFAULT_SHORTEST_PATH_THRESHOLD = 1;
    private static final double DEFAULT_PVALUE_CUT_OFF = 0.05;
    private static final double DEFAULT_MODULE_FREQUENCY_THRESHOLD = 0.1;
    private static final double DEFAULT_MRNA_ZSCORE_THRESHOLD = 1.96;
    private static final int DEFAULT_SEED_RADIUS_THRESHOLD = 2;
    private static final int DEFAULT_NUM_RANDOM_TRIALS = 1000;
    private static final double DEFAULT_GENE_FREQUENCY_THRESHOLD = 0.00001;

    private boolean veryQuiet;
    private int numLocalTrials;
    private int numGlobalTrials;
    private File backgroundGeneFile;
    private File geneFile;
    private boolean debugMode = false;
    private int shortestPathThreshold;
    private double pValueCutOff;
    private boolean interactiveMode = false;
    private File mutationFile;
    private File cnaFile;
    private File mRNAFile;
    private File caseFile;
    private double deltaThreshold;
    private double moduleFrequencyThreshold;
    private double mRNAZScoreThreshold;
    private double geneFrequencyThreshold;
    private int seedRadiusThreshold;
    private String title;
    private boolean includeLowLevelCnaChanges = false;
    private String networkPartitionAlgorithm;
    private boolean identifyModules = true;

    private static final String GENE_FILE = "gene_file";
    private static final String MUTATION_FILE = "mutation_file";
    private static final String CNA_FILE = "cna_file";
    private static final String MRNA_FILE = "mrna_file";
    private static final String CASE_FILE = "case_file";
    private static final String TITLE = "title";
    private static final String MODULE_FREQUENCY_THRESHOLD = "module_frequency_threshold";
    private static final String DELTA_THRESHOLD = "delta_threshold";
    private static final String SEED_RADIUS_THRESHOLD = "seed_radius_threshold";
    private static final String SHORTEST_PATH_THRESHOLD = "shortest_path_threshold";
    private static final String NUM_LOCAL_TRIALS = "num_local_trials";
    private static final String NUM_GLOBAL_TRIALS = "num_global_trials";
    private static final String MRNA_ZSCORE_THRESHOLD = "mrna_zscore_threshold";
    private static final String GENE_FREQUENCY_THRESHOLD = "gene_frequency_threshold";
    private static final String NETWORK_PARTITION_ALGORITHM = "network_partition_algorithm";
    private static final String P_VALUE_THRESHOLD = "p_value_threshold";
    private static final String IDENTIFY_MODULES = "identify_modules";
    private static final String INCLUDE_LOW_LEVEL_CNA_CHANGES
            = "include_low_level_cna_changes";

    /**
     * Newman-Girvan Algorithm.
     */
    public static final String NG = "ng";

    /**
     * Simulated Annealing Algorithm.
     */
    public static final String SA = "sa";

    /**
     * Private Constructor.  Enforces Singleton Pattern.
     */
    private GlobalConfig() {
        this.numLocalTrials = DEFAULT_NUM_RANDOM_TRIALS;
        this.shortestPathThreshold = DEFAULT_SHORTEST_PATH_THRESHOLD;
        this.pValueCutOff = DEFAULT_PVALUE_CUT_OFF;
        this.moduleFrequencyThreshold = DEFAULT_MODULE_FREQUENCY_THRESHOLD;
        this.seedRadiusThreshold = DEFAULT_SEED_RADIUS_THRESHOLD;
        this.mRNAZScoreThreshold = DEFAULT_MRNA_ZSCORE_THRESHOLD;
        this.geneFrequencyThreshold = DEFAULT_GENE_FREQUENCY_THRESHOLD;
        this.networkPartitionAlgorithm = NG;
        this.numLocalTrials = 0;
        this.numGlobalTrials = 0;
    }

    /**
     * Gets Singleton Instance.
     *
     * @return Singleton GlobalConfig Instance.
     */
    public static GlobalConfig getInstance() {
        if (globalConfig == null) {
            globalConfig = new GlobalConfig();
        }
        return globalConfig;
    }

    /**
     * Loads Properties from the Specified File.
     *
     * @param file Properties File.
     * @throws IOException IO Error.
     */
    public void loadProperties(File file) throws IOException {
        Properties properties = new Properties();
        properties.load(new FileInputStream(file));

        String value = properties.getProperty(GENE_FILE);
        if (value != null) {
            this.geneFile = new File(value);
        }

        value = properties.getProperty(MUTATION_FILE);
        if (value != null) {
            this.mutationFile = new File(value);
        }

        value = properties.getProperty(CNA_FILE);
        if (value != null) {
            this.cnaFile = new File(value);
        }

        value = properties.getProperty(MRNA_FILE);
        if (value != null) {
            this.mRNAFile = new File(value);
        }

        value = properties.getProperty(CASE_FILE);
        if (value != null) {
            this.caseFile = new File(value);
        }

        this.title = properties.getProperty(TITLE);

        value = properties.getProperty(MODULE_FREQUENCY_THRESHOLD);
        if (value != null) {
            try {
                this.moduleFrequencyThreshold = Double.parseDouble(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  " + MODULE_FREQUENCY_THRESHOLD
                        + " must be a number.");
            }
        }

        value = properties.getProperty(DELTA_THRESHOLD);
        if (value != null) {
            try {
                this.deltaThreshold = Double.parseDouble(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  " + DELTA_THRESHOLD
                        + " must be a number.");
            }
        }

        value = properties.getProperty(MRNA_ZSCORE_THRESHOLD);
        if (value != null) {
            try {
                this.mRNAZScoreThreshold = Double.parseDouble(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  " + MRNA_ZSCORE_THRESHOLD
                        + " must be a number.");
            }
        }

        value = properties.getProperty(GENE_FREQUENCY_THRESHOLD);
        if (value != null) {
            try {
                this.geneFrequencyThreshold = Double.parseDouble(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  " + GENE_FREQUENCY_THRESHOLD
                        + " must be a number.");
            }
        }

        value = properties.getProperty(SEED_RADIUS_THRESHOLD);
        if (value != null) {
            try {
                this.seedRadiusThreshold = Integer.parseInt(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  " + SEED_RADIUS_THRESHOLD
                        + " must be an interger.");
            }
        }

        value = properties.getProperty(SHORTEST_PATH_THRESHOLD);
        if (value != null) {
            try {
                this.shortestPathThreshold = Integer.parseInt(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  " + SHORTEST_PATH_THRESHOLD
                        + " must be an interger.");
            }
        }

        value = properties.getProperty(P_VALUE_THRESHOLD);
        if (value != null) {
            try {
                this.pValueCutOff = Double.parseDouble(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  " + P_VALUE_THRESHOLD
                        + " must be a number.");
            }
        }

        value = properties.getProperty(NUM_LOCAL_TRIALS);
        if (value != null) {
            try {
                this.numLocalTrials = Integer.parseInt(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  " + NUM_LOCAL_TRIALS
                        + " must be an interger.");
            }
        }

        value = properties.getProperty(NUM_GLOBAL_TRIALS);
        if (value != null) {
            try {
                this.numGlobalTrials = Integer.parseInt(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  " + NUM_GLOBAL_TRIALS
                        + " must be an interger.");
            }
        }

        value = properties.getProperty(INCLUDE_LOW_LEVEL_CNA_CHANGES);
        if (value != null) {
            try {
                this.includeLowLevelCnaChanges = Boolean.parseBoolean(value);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Property:  "
                        + INCLUDE_LOW_LEVEL_CNA_CHANGES
                        + " must be set to true or false.");
            }
        }

        value = properties.getProperty(NETWORK_PARTITION_ALGORITHM);
        if (value != null) {
            if (value.equals(NG) || value.equals(SA)) {
                this.networkPartitionAlgorithm = value;
            } else {
                throw new IllegalArgumentException("Unsupported option:  " + value + " for "
                        + "property:  " + NETWORK_PARTITION_ALGORITHM + ".");
            }
        }

        value = properties.getProperty(IDENTIFY_MODULES);
        if (value != null) {
            if (value.equals("F") || value.equals ("FALSE")) {
                identifyModules = false;
            }
        }

        if (geneFile == null) {
            if (caseFile == null) {
                throw new IllegalArgumentException("Properties file is missing:  " + CASE_FILE
                        + " property.");
            }
        }
    }

    /**
     * Gets the Very Quiet Flag.
     *
     * @return very quiet flag.
     */
    public boolean isVeryQuiet() {
        return veryQuiet;
    }

    /**
     * Sets the Very Quiet Flag.
     *
     * @param q Very Quiet Flag.
     */
    public void setBeVeryQuiet(boolean q) {
        this.veryQuiet = q;
    }

    /**
     * Gets Number of Random Trials.
     *
     * @return number of random trials.
     */
    public int getNumLocalTrials() {
        return numLocalTrials;
    }

    /**
     * Sets the Number of Random Trials.
     *
     * @param n Number of Random Trials.
     */
    public void setNumLocalTrials(int n) {
        this.numLocalTrials = n;
    }

    /**
     * Gets the Background Gene File.
     *
     * @return Background Gene File.
     */
    public File getBackgroundGeneFile() {
        return backgroundGeneFile;
    }

    /**
     * Sets the Background Gene File.
     *
     * @param f Background Gene File.
     */
    public void setBackgroundGeneFile(File f) {
        this.backgroundGeneFile = f;
    }

    /**
     * Returns true if running in debug mode.
     *
     * @return debug mode.
     */
    public boolean isDebugMode() {
        return debugMode;
    }

    /**
     * Sets the debug mode.
     *
     * @param d debug mode.
     */
    public void setDebugMode(boolean d) {
        this.debugMode = d;
    }

    /**
     * Gets the Shortest Path Threshold.
     *
     * @return gets the shortest path threshold.
     */
    public int getShortestPathThreshold() {
        return shortestPathThreshold;
    }

    /**
     * Sets the Shortest Path Threshold.
     *
     * @param s Shortest Path Threshold.
     */
    public void setShortestPathThreshold(int s) {
        this.shortestPathThreshold = s;
    }

    /**
     * Gets the P-Value Cutoff.
     *
     * @return p-value cutoff.
     */
    public double getPValueCutOff() {
        return pValueCutOff;
    }

    /**
     * Sets the P-Value Cutoff.
     *
     * @param p p-value cutoff.
     */
    public void setPValueCutOff(double p) {
        this.pValueCutOff = p;
    }

    /**
     * Gets the Default Shortest Path Threshold.
     *
     * @return default shorest path threshold.
     */
    public int getDefaultShortestPathThreshold() {
        return DEFAULT_SHORTEST_PATH_THRESHOLD;
    }

    /**
     * Gets the Default PValue CutOff.
     *
     * @return default pvalue cutoff.
     */
    public double getDefaultPValueCutOff() {
        return DEFAULT_PVALUE_CUT_OFF;
    }

    /**
     * Gets the Gene File.
     *
     * @return Gene File.
     */
    public File getGeneFile() {
        return geneFile;
    }

    /**
     * Sets the Gene File.
     *
     * @param f Gene File.
     */
    public void setGeneFile(File f) {
        this.geneFile = f;
    }

    /**
     * Currently running in interactive mode.
     *
     * @return returns true if running in interactive mode.
     */
    public boolean isInteractiveMode() {
        return interactiveMode;
    }

    /**
     * Sets the interactive mode flag.
     *
     * @param m interactive mode flag.
     */
    public void setInteractiveMode(boolean m) {
        this.interactiveMode = m;
    }

    /**
     * Gets the Mutation File.
     *
     * @return Mutation File.
     */
    public File getMutationFile() {
        return mutationFile;
    }

    /**
     * Gets the CNA File.
     *
     * @return CNA File.
     */
    public File getCnaFile() {
        return cnaFile;
    }

    /**
     * Gets the mRNA File.
     *
     * @return mRNA File.
     */
    public File getMRNAFile() {
        return mRNAFile;
    }

    /**
     * Gets the Case File.
     *
     * @return Case File.
     */
    public File getCaseFile() {
        return caseFile;
    }

    /**
     * Gets the Module Frequency Threshold.
     *
     * @return Module Frequency Threshold.
     */
    public double getModuleFrequencyThreshold() {
        return moduleFrequencyThreshold;
    }

    /**
     * Gets the Delta Threshold for Adding New Genes to a Module.
     *
     * @return Delta Threshold.
     */
    public double getDeltaThreshold() {
        return deltaThreshold;
    }

    /**
     * Gets the mRNA Z-Score Threshold.
     *
     * @return mRNA Z-Score Threshold.
     */
    public double getMRNAZScoreThreshold() {
        return mRNAZScoreThreshold;
    }

    /**
     * Gets the Seed Radius Threshold.
     *
     * @return seed radius threshold.
     */
    public int getSeedRadiusThreshold() {
        return seedRadiusThreshold;
    }

    /**
     * Gets the Gene Frequency Threshold.
     *
     * @return Gene Frequency Threshold.
     */
    public double getGeneFrequencyThreshold() {
        return geneFrequencyThreshold;
    }

    /**
     * Gets Report Title.
     *
     * @return Report Title.
     */
    public String getTitle() {
        return title;
    }

    /**
     * Include Hemizygous Deletions and Single Copy Gains.
     *
     * @return true of false.
     */
    public boolean includeLowLevelCnaChanges() {
        return includeLowLevelCnaChanges;
    }

    /**
     * Gets the user-specified Network Partioning Algorithm.
     *
     * @return network partitioning algorithm.
     */
    public String getNetworkPartitionAlgorithm() {
        return networkPartitionAlgorithm;
    }

    /**
     * Sets the user-specified Network Partitioning Algorithm.
     *
     * @param algo Network partitioning algorithm.
     */
    public void setNetworkPartitionAlgorithm(String algo) {
        this.networkPartitionAlgorithm = algo;
    }

    /**
     * Gets the number of trials to execute in the Global Random Null Model.
     *
     * @return number of trials.
     */
    public int getNumGlobalTrials() {
        return numGlobalTrials;
    }

    /**
     * Sets the number of trials to execute in the Global Random Null Model.
     *
     * @param n number of trials.
     */
    public void setNumGlobalTrials(int n) {
        this.numGlobalTrials = n;
    }

    /**
     * Gets whether to identify modules.
     * @return
     */
    public boolean identifyModules() {
        return identifyModules;
    }
}
