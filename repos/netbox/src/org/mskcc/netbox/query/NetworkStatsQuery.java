package org.mskcc.netbox.query;

import org.hibernate.Query;
import org.hibernate.Session;
import org.mskcc.netbox.model.NetworkStats;
import org.mskcc.netbox.util.GlobalSession;

/**
 * Network Stats Queries.
 *
 * @author Ethan Cerami.
 */
public final class NetworkStatsQuery {

    /**
     * Private Constructor to prevent instantiation.
     */
    private NetworkStatsQuery() {
    }

    /**
     * Gets the Current Network Stats from the Database.
     *
     * @return NetworkStats Object.
     */
    public static NetworkStats getNetworkStats() {
        Session session = GlobalSession.getInstance().getSession();
        Query query = session.getNamedQuery("org.mskcc.netbox.getNetworkStats");
        return (NetworkStats) query.uniqueResult();
    }

    /**
     * Deletes all Network Stats from the Database.
     */
    public static void deleteAllNetworkStats() {
        Session session = GlobalSession.getInstance().getSession();
        Query query = session.getNamedQuery("org.mskcc.netbox.deleteAllNetworkStats");
        query.executeUpdate();
    }
}
