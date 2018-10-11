package org.mskcc.netbox.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * Misc. File Utilities.
 *
 * @author Ethan Cerami.
 */
public final class FileUtil {

    /**
     * Private constructor, to prevent instantiation.
     */
    private FileUtil() {
    }

    /**
     * Gets Number of Lines in Specified File.
     *
     * @param file File.
     * @return number of lines.
     * @throws java.io.IOException Error Reading File.
     */
    public static int getNumLines(File file) throws IOException {
        int numLines = 0;
        FileReader reader = new FileReader(file);
        BufferedReader buffered = new BufferedReader(reader);
        String line = buffered.readLine();
        while (line != null) {
            if (!line.startsWith("#") && line.trim().length() > 0) {
                numLines++;
            }
            line = buffered.readLine();
        }
        reader.close();
        return numLines;
    }

    /**
     * Gets Next Line of Input.  Filters out Empty Lines and Comments.
     *
     * @param buf BufferedReader Object.
     * @return next line of input.
     * @throws IOException Error reading input stream.
     */
    public static String getNextLine(BufferedReader buf) throws IOException {
        String line = buf.readLine();
        while (line != null && (line.trim().length() == 0
                || line.trim().startsWith("#"))) {
            line = buf.readLine();
        }
        return line;
    }

    /**
     * Tests if the Specified File Exists.  Aborts command line if file does not exist.
     *
     * @param file Target File.
     */
    public static void testFileExists(File file) {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(file));
        } catch (FileNotFoundException e) {
            CommandLineUtil.abort("File Not Found:  " + file.getAbsolutePath());
        } catch (IOException e) {
            CommandLineUtil.abort("Error Occurred while trying to read file:  "
                    + file.getAbsolutePath());
        } finally {
            try {
                reader.close();
            } catch (IOException e) {
                reader = null;
            }
        }
    }
}
