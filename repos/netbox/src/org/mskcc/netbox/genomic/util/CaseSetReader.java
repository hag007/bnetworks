package org.mskcc.netbox.genomic.util;

import org.mskcc.netbox.util.CommandLineUtil;
import org.mskcc.netbox.util.FileUtil;
import org.mskcc.netbox.util.ProgressMonitor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

/**
 * Utility Class for Reading in a Set of Case IDs from a specified file.
 *
 * @author Ethan Cerami.
 */
public final class CaseSetReader {
    private HashSet<String> caseIdSet = new HashSet<String>();

    /**
     * Constructor.
     *
     * @param file File containing case IDs.
     * @throws IOException IO Error.
     */
    public CaseSetReader(File file) throws IOException {
        ProgressMonitor pMonitor = ProgressMonitor.getInstance();
        pMonitor.setCurrentMessage("Reading cases from:  " + file.getAbsolutePath());
        pMonitor.setMaxValue(FileUtil.getNumLines(file));
        FileReader reader = new FileReader(file);
        BufferedReader buf = new BufferedReader(reader);

        String line = buf.readLine();
        while (line != null) {
            line = line.trim();
            if (!line.startsWith("#")) {
                caseIdSet.add(line);
            }
            line = buf.readLine();
            pMonitor.incrementCurValue();
            CommandLineUtil.showProgress(pMonitor);
        }
    }

    /**
     * Gets the Set of Case IDs.
     *
     * @return Hashset of Case IDs.
     */
    public HashSet<String> getCaseIdSet() {
        return caseIdSet;
    }
}
