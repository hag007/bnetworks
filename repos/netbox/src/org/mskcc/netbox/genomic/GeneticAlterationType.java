package org.mskcc.netbox.genomic;

/**
 * Genetic Aleration Type Enumerator.
 */
public final class GeneticAlterationType {
    private String type;
    /**
     * Homozygous Deletion.
     */
    public static final String HOMOZYGOUS_DELETION = "-2";

    /**
     * Hemizygous Deletion.
     */
    public static final String HEMIZYGOUS_DELETION = "-1";

    /**
     * Zero Change.
     */
    public static final String ZERO = "0";

    /**
     * Single Copy Gain.
     */
    public static final String GAIN = "1";

    /**
     * Multi-Copy Amplification.
     */
    public static final String AMPLIFICATION = "2";

    /**
     * Data Not Available.
     */
    public static final String NAN = "NaN";

    /**
     * Private Constructor. Enumeration Pattern.
     *
     * @param t Alteration Type.
     */
    private GeneticAlterationType(String t) {
        this.type = t;
    }

    /**
     * Gets Type Name.
     *
     * @return Type Name.
     */
    public String toString() {
        return type;
    }

    /**
     * Get Type by Type Name.
     *
     * @param type Type Name, e.g. "MUTATION" or "COPY_NUMBER_ALERATION" or "MRNA_EXPRESSION"".
     * @return correct GeneticAlterationType Object.
     */
    public static GeneticAlterationType getType(String type) {
        if (type == null) {
            throw new IllegalArgumentException("type is null");
        } else if (type.equals(MUTATION.toString())) {
            return MUTATION;
        } else if (type.equals(COPY_NUMBER_ALTERATION.toString())) {
            return COPY_NUMBER_ALTERATION;
        } else if (type.equals(MRNA_EXPRESSION.toString())) {
            return MRNA_EXPRESSION;
        } else {
            throw new NullPointerException("Cannot find:  " + type);
        }
    }

    /**
     * Type:  MUTATION.
     */
    public static final GeneticAlterationType MUTATION
            = new GeneticAlterationType("MUTATION");

    /**
     * Type:  COPY_NUMBER_ALTERATION.
     */
    public static final GeneticAlterationType COPY_NUMBER_ALTERATION
            = new GeneticAlterationType("COPY_NUMBER_ALTERATION");

    /**
     * Type:  MRNA_EXPRESSION.
     */
    public static final GeneticAlterationType MRNA_EXPRESSION
            = new GeneticAlterationType("MRNA_EXPRESSION");

    /**
     * Type:  MERGED_PROFILE.
     */
    public static final GeneticAlterationType MERGED_PROFILE
            = new GeneticAlterationType("MERGED_PROFILE");

}

