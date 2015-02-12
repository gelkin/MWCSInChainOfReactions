package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import ec.EvolutionState;

import java.awt.font.GraphicAttribute;
import java.io.*;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class Solver {
    private static final String GENOME_SIZE_FROM_FILE = "genome-size-from-file";
    private static final String FILE_PARAM = "-file";
    private static final String EVOLVE_SUFFIX = "-evolve.";
    private static final String STAT_SUFFIX = ".stat";
    private static final String VERTICES_OUTFILE_NAME = "./src/main/resources/nodes.tsv";
    private static final String EDGES_OUTFILE_NAME = "./src/main/resources/edges.tsv";
    private static final int STRINGS_TO_SKIP = 5;

    private Graph graph;

    public static void main(String[] args) {
        (new Solver()).run(args);
    }

    private void run(String[] args) {
        if (args.length < 4) {
            System.out.println("Too few arguments");
            return;
        }

        String verticesFile = args[0];
        String edgesFile = args[1];
        String signalsFile = args[2];
        String parametersFile = args[3];

        graph = initGraphFromFiles(verticesFile, edgesFile, signalsFile);
        MWCGProblem.setGraph(graph);

        String newParametersFile = "";
        try {
            int numberOfEdges = getNumberOfLines(edgesFile);

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

        final String[] newArgs = new String[2];
        newArgs[0] = FILE_PARAM;
        newArgs[1] = newParametersFile;

        // If you want get .tsv files from mwcs-single.stat, just comment next line
        // ( otherwise ec.Evolve.main(...) calls System.exit(0) )
        // ec.Evolve.main(newArgs);

        // Let's return result subgraphin readable form
        String[] parts = parametersFile.split("\\.");
        String statFile = parts[0] + STAT_SUFFIX;
        List<Boolean> bestIndividual = getBestIndividual(statFile);

        int heaviestComponentNumber = (int) graph.getHeaviestComponentInfo(bestIndividual)[1];
        List<Boolean> heaviestComponent = graph.getComponentByNumber(bestIndividual, heaviestComponentNumber);
        writeResults(heaviestComponent);
    }

    private Graph initGraphFromFiles(String verticesFile, String edgesFile, String signalsFile) {

        // Vertices:
        LinkedHashMap<String, String> verticesToSignal = new LinkedHashMap<>();
        try(BufferedReader br = new BufferedReader(new FileReader(verticesFile))) {
            String line = br.readLine();

            while (line != null) {
                String[] vertexInfo = line.split("\\t");
                String vertex = vertexInfo[0];
                String signal = vertexInfo[1];

                verticesToSignal.put(vertex, signal);

                line = br.readLine();
            }

        } catch (FileNotFoundException e) {
            System.out.println("Cannot find file " + verticesFile + ".");
            e.printStackTrace();
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file " + verticesFile + ".");
            e.printStackTrace();
        }

        // Edges:
        LinkedHashMap<String, LinkedHashMap<String, String>> edgesAsMap = new LinkedHashMap<>();
        List<Edge<String>> edges = new ArrayList<>();
        Map<String, Map<String, Integer>> edgesToIndex = new HashMap<>();

        int edgeNumber = 0;
        try(BufferedReader br = new BufferedReader(new FileReader(edgesFile))) {
            String line = br.readLine();

            while (line != null) {
                String[] edgeInfo = line.split("\\t");
                String from = edgeInfo[0];
                String to = edgeInfo[1];
                String signal = edgeInfo[2];

                if (!edgesAsMap.containsKey(from)) {
                    edgesAsMap.put(from, new LinkedHashMap<String, String>());
                    edgesToIndex.put(from, new HashMap<String, Integer>());
                }

                if (!edgesAsMap.containsKey(to)) {
                    edgesAsMap.put(to, new LinkedHashMap<String, String>());
                    edgesToIndex.put(to, new HashMap<String, Integer>());
                }


                edgesAsMap.get(from).put(to, signal);
                edgesAsMap.get(to).put(from, signal);

                edges.add(new Edge<>(from, to));

                edgesToIndex.get(from).put(to, edgeNumber);
                edgesToIndex.get(to).put(from, edgeNumber);

                ++edgeNumber;

                line = br.readLine();
            }

        } catch (FileNotFoundException e) {
            System.out.println("Cannot find file " + edgesFile + ".");
            e.printStackTrace();
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file " + edgesFile + ".");
            e.printStackTrace();
        }

        // Signals:
        Map<String, Double> signals = new HashMap<>();
        try(BufferedReader br = new BufferedReader(new FileReader(signalsFile))) {
            String line = br.readLine();

            while (line != null) {
                String[] vertexInfo = line.split("\\t");
                String signal = vertexInfo[0];
                double value = Double.parseDouble(vertexInfo[1]);

                signals.put(signal, value);

                line = br.readLine();
            }

        } catch (FileNotFoundException e) {
            System.out.println("Cannot find file " + verticesFile + ".");
            e.printStackTrace();
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file " + verticesFile + ".");
            e.printStackTrace();
        }

        return new Graph<>(edgesAsMap, verticesToSignal, signals, edges, edgesToIndex);
    }

    private int getNumberOfLines(String fileName) throws IOException {
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
    private String createFileWithReplacedString(String fileName, String oldString, String newString) throws IOException {
        Path oldPath = Paths.get(fileName);
        Charset charset = StandardCharsets.UTF_8;

        String content = new String(Files.readAllBytes(oldPath), charset);
        content = content.replaceAll(oldString, newString);

        String[] parts = fileName.split("\\.");
        String newFile = parts[0] + EVOLVE_SUFFIX + parts[1];
        Path newPath = Paths.get(newFile);
        Files.write(newPath, content.getBytes(charset));
        return newFile;
    }

    private List<Boolean> getBestIndividual(String fileName) {
        String individual = "";
        try (FileInputStream fs = new FileInputStream(fileName)) {
            BufferedReader br = new BufferedReader(new InputStreamReader(fs));
            for(int i = 0; i < STRINGS_TO_SKIP; ++i) {
                br.readLine();
            }
            individual = br.readLine();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        if (individual != "") {
            List<Boolean> mask = new ArrayList<>(individual.length());
            for (int i = 0; i < individual.length(); ++i) {
                if (individual.charAt(i) == '0') {
                    mask.add(false);
                } else {
                    mask.add(true);
                }
            }

            return mask;
        }

        return null;
    }

    private void writeResults(List<Boolean> mask) {
        // vertices
        List<Pair<String, String>> verticesToSignals = graph.getVerticesToSignalsByMask(mask);
        try (FileOutputStream fos = new FileOutputStream(VERTICES_OUTFILE_NAME)) {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));

            for (int i = 0; i < verticesToSignals.size(); i++) {
                Pair<String, String> vertex = verticesToSignals.get(i);
                bw.write(vertex.first + "\t" + vertex.second + "\t" + "Node " + i);
                bw.newLine();
            }

            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // edges
        List<Pair<Edge<String>, String>> edgesToSignals = graph.getEdgesToSignalsByMask(mask);
        try (FileOutputStream fos = new FileOutputStream(EDGES_OUTFILE_NAME)) {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));

            for (int i = 0; i < edgesToSignals.size(); i++) {
                Pair<Edge<String>, String> edge = edgesToSignals.get(i);
                bw.write(i + "\t" + edge.first.first + "\t" + edge.first.second + "\t" + edge.second + "\t" + "Edge " + i);
                bw.newLine();
            }

            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
