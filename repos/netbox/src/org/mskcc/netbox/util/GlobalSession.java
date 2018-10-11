package org.mskcc.netbox.util;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.AnnotationConfiguration;
import org.hibernate.cfg.Configuration;

import java.io.File;

/**
 * Global Session Class.
 *
 * @author Ethan Cerami
 */
public final class GlobalSession {
    private static GlobalSession globalSession;
    private SessionFactory sessionFactory;
    private Session session;

    /**
     * Private Constructor to Enforce Singleton Pattern.
     */
    private GlobalSession() {
        String netBoxHome = Constants.getNetBoxHome();
        File conf = new File(netBoxHome + "/config/hibernate.cfg.xml");
        Configuration config = new AnnotationConfiguration().configure(conf);
        config.setProperty("hibernate.connection.url", "jdbc:hsqldb:" + netBoxHome + "/db/netbox");
        sessionFactory = config.buildSessionFactory();
        session = sessionFactory.openSession();
    }

    /**
     * Gets Singleton Instance.
     *
     * @return Global Session Object.
     */
    public static GlobalSession getInstance() {
        if (globalSession == null) {
            globalSession = new GlobalSession();
        }
        return globalSession;
    }

    /**
     * Gets the Session Factory.
     *
     * @return Session Factory.
     */
    public SessionFactory getSessionFactory() {
        return sessionFactory;
    }

    /**
     * Gets the Session.
     *
     * @return Session Object.
     */
    public Session getSession() {
        return session;
    }

    /**
     * Closes Session Factory and Session.
     */
    public void closeAll() {
        //  Must close session before closing session factory.
        session.close();
        sessionFactory.close();
    }
}
