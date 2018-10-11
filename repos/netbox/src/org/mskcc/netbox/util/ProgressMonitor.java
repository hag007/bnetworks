package org.mskcc.netbox.util;

import java.util.ArrayList;

/**
 * Monitors Progress of Long Term Tasks.
 *
 * @author Ethan Cerami.
 */
public final class ProgressMonitor {
    private int maxValue;
    private int curValue;
    private String currentMessage;
    private StringBuffer log = new StringBuffer();
    private boolean consoleMode;
    private ArrayList<String> warningList = new ArrayList<String>();
    private static ProgressMonitor pMonitor;

    private ProgressMonitor() {
    }

    /**
     * Gets the Progress Monitor Singleton.
     *
     * @return Progress Monitor Singleton.
     */
    public static ProgressMonitor getInstance() {
        if (pMonitor == null) {
            pMonitor = new ProgressMonitor();
        }
        return pMonitor;
    }

    /**
     * Sets Console Flag.
     * When set to true Progress Monitor Messages are displayed to System.out.
     *
     * @param consoleFlag Console Mode Flag.
     */
    public void setConsoleMode(boolean consoleFlag) {
        this.consoleMode = consoleFlag;
    }

    /**
     * Gets Console Mode Flag.
     *
     * @return Boolean Flag.
     */
    public boolean isConsoleMode() {
        return this.consoleMode;
    }

    /**
     * Gets Percentage Complete.
     *
     * @return double value.
     */
    public double getPercentComplete() {
        if (curValue == 0) {
            return 0.0;
        } else {
            return (curValue / (double) maxValue);
        }
    }

    /**
     * Gets Max Value.
     *
     * @return max value.
     */
    public int getMaxValue() {
        return maxValue;
    }

    /**
     * Sets Max Value.
     *
     * @param m Max Value.
     */
    public void setMaxValue(int m) {
        this.maxValue = m;
        this.curValue = 0;
    }

    /**
     * Gets Current Value.
     *
     * @return Current Value.
     */
    public int getCurValue() {
        return curValue;
    }

    /**
     * Incremenets the Current Value.
     */
    public void incrementCurValue() {
        curValue++;
    }

    /**
     * Sets the Current Value.
     *
     * @param value Current Value.
     */
    public void setCurValue(int value) {
        this.curValue = value;
    }

    /**
     * Gets the Current Task Message.
     *
     * @return Currest Task Message.
     */
    public String getCurrentMessage() {
        return currentMessage;
    }

    /**
     * Gets Log of All Messages.
     *
     * @return String Object.
     */
    public String getLog() {
        return log.toString();
    }

    /**
     * Logs a Message.
     *
     * @param s Current Task Message.
     */
    public void setCurrentMessage(String s) {
        this.currentMessage = s;
        this.log.append(currentMessage + "\n");
        if (consoleMode && !GlobalConfig.getInstance().isVeryQuiet()) {
            System.err.println(currentMessage);
        }
    }

    /**
     * Logs a warning message.
     *
     * @param warning warning message.
     */
    public void logWarning(String warning) {
        warningList.add(warning);
    }

    /**
     * Gets a list of warning messages.
     *
     * @return ArrayList of warning messages.
     */
    public ArrayList<String> getWarnings() {
        return warningList;
    }
}
