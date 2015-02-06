package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class Solver {
    private static final String GENOME_SIZE_FROM_FILE = "genome-size-from-file";
    private static final String FILE_PARAM = "-file";

    public static void main(String[] args) {
        if (args.length < 4) {
            System.out.println("Too few arguments");
            return;
        }

        String verticesFile = args[0];
        String edgesFile = args[1];
        String signalsFile = args[2];
        String parametersFile = args[3];

        MWCGProblem.setSourceFiles(verticesFile, edgesFile, signalsFile);

        String newParametersFile = "";
        try {
            int numberOfEdges = getNumberOfLines(edgesFile);
            // todo
            System.out.println("numberOfEdges = " + numberOfEdges);

            if (numberOfEdges == -1) {
                return;
            }

            newParametersFile = createFileWithReplacedString(parametersFile,
                                                             GENOME_SIZE_FROM_FILE,
                                                             (new Integer(numberOfEdges)).toString());
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file.");
            e.printStackTrace();
        }

        String[] newArgs = new String[2];
        newArgs[0] = FILE_PARAM;
        newArgs[1] = newParametersFile;

        ec.Evolve.main(newArgs);
    }

    private static int getNumberOfLines(String fileName) throws IOException {
        int linesNumber = 0;
        try(BufferedReader br = new BufferedReader(new FileReader(fileName))) {
            String line = br.readLine();

            while (line != null) {
                ++linesNumber;
                line = br.readLine();
            }
        } catch (FileNotFoundException e) {
            System.out.println("Cannot find file " + fileName + ".");
            linesNumber = -1;
            e.printStackTrace();
        }

        return linesNumber;
    }

    // fileName.ext -> fileName-evolve.ext
    private static String createFileWithReplacedString(String fileName, String oldString, String newString) throws IOException {
        Path oldPath = Paths.get(fileName);
        Charset charset = StandardCharsets.UTF_8;

        String content = new String(Files.readAllBytes(oldPath), charset);
        content = content.replaceAll(oldString, newString);

        String[] parts = fileName.split("\\.");
        String newFile = parts[0] + "-evolve." + parts[1];
        Path newPath = Paths.get(newFile);
        Files.write(newPath, content.getBytes(charset));
        return newFile;
    }
}
