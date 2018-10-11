package org.mskcc.netbox.util;

/**
 * Global Constants.
 *
 * @author Ethan Cerami.
 */
public final class Constants {

    /**
     * Private Constructor to prevent instantiation.
     */
    private Constants() {
    }

    /**
     * Gets NETBOX_HOME Environment Variable.
     *
     * @return NetBox Home.
     */
    public static String getNetBoxHome() {
        String netBoxHome = System.getenv("NETBOX_HOME");
        if (netBoxHome == null) {
            System.out.println("$NETBOX_HOME is not set.  Aborting.");
            System.exit(0);
        }
        return netBoxHome;
    }
}
