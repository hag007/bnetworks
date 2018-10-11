package org.mskcc.netbox.util;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;

/**
 * Command Line Utilities.
 *
 * @author Ethan Cerami
 */
public final class CommandLineUtil {

    /**
     * Private Constructor to prevent instantiation.
     */
    private CommandLineUtil() {
    }

    /**
     * Aborts the Command Line Program with the User Error Message.
     *
     * @param userMsg User Message.
     */
    public static void abort(String userMsg) {
        System.out.println(userMsg);
        System.exit(1);
    }

    /**
     * Outputs Progress to Console.
     *
     * @param pMonitor ProgressMonitor Object.
     */
    public static synchronized void showProgress(ProgressMonitor pMonitor) {
        String msg;
        if (pMonitor.isConsoleMode()) {
            int currentValue = pMonitor.getCurValue();
            System.err.print(".");
            if (currentValue % 100 == 0) {
                NumberFormat format = DecimalFormat.getPercentInstance();
                double percent = pMonitor.getPercentComplete();
                msg = new String("Percentage Complete:  "
                        + format.format(percent) + ":  " + pMonitor.getCurrentMessage());
                System.err.println("\n" + msg);
                Runtime rt = Runtime.getRuntime();
                long used = rt.totalMemory() - rt.freeMemory();
                System.err.println("Mem Allocated:  " + getMegabytes(rt.totalMemory())
                        + ", Mem used:  " + getMegabytes(used) + ", Mem free:  "
                        + getMegabytes(rt.freeMemory()));
            }
            if (currentValue == pMonitor.getMaxValue()) {
                System.err.println();
            }
        }
    }

    /**
     * Shows Warning Messages.
     *
     * @param pMonitor Progress Monitor Object.
     */
    public static void showWarnings(ProgressMonitor pMonitor) {
        ArrayList warningList = pMonitor.getWarnings();
        if (warningList.size() == 0) {
            System.err.println("\nNo warning/error messages generated.");
        } else {
            System.err.println("\nWarnings / Errors:");
            System.err.println("-------------------");
            for (int i = 0; i < warningList.size(); i++) {
                System.err.println(i + ".  " + warningList.get(i));
            }
        }
    }

    private static String getMegabytes(long bytes) {
        double mBytes = (bytes / 1024.0) / 1024.0;
        DecimalFormat formatter = new DecimalFormat("#,###,###.###");
        return formatter.format(mBytes) + " MB";
    }
}
