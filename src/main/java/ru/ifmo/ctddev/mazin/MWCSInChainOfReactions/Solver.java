package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import java.io.*;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Solver {
    private static final String FILE_PARAM = "-file";
    private static final String PARAM_IDENTIFIER = "-p";
    private static final String GENOME_SIZE_PARAM = "pop.subpop.0.species.genome-size=";
    private static final String VERTICES_OUTFILE_NAME = "./src/main/resources/nodes.tsv";
    private static final String EDGES_OUTFILE_NAME = "./src/main/resources/edges.tsv";
    private static final String NUMBER_TO_FITNESS_CSV = "./src/main/resources/number_to_fitness.csv";
    private String SOLUTION_COMPONENTS_FILE = "./src/main/resources/solution-components.stat";

    private String parametersFile = "./src/main/resources/mwcs-multi.params";
    private String statFile = "./src/main/resources/mwcs-multi.stat";

    private Graph graph;

    public static void main(String[] args) {
        (new Solver()).run(args);
    }

    private void run(String[] args) {
        if (args.length < 3) {
            System.out.println("Too few arguments");
            return;
        }

        String verticesFile = args[0];
        String edgesFile = args[1];
        String signalsFile = args[2];

        graph = initGraphFromFiles(verticesFile, edgesFile, signalsFile);
        MWCGProblem.setGraph(graph);
        GraphMutatorPipeline.setGraph(graph);
        GraphIndividual.setGraph(graph);

        int numberOfEdges;
        try {

            numberOfEdges = getNumberOfLines(edgesFile);
            if (numberOfEdges == -1) {
                return;
            }
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file.");
            e.printStackTrace();
            return;
        }

        final String[] newArgs = new String[4];
        newArgs[0] = FILE_PARAM;
        newArgs[1] = parametersFile;
        newArgs[2] = PARAM_IDENTIFIER;
        newArgs[3] = GENOME_SIZE_PARAM + numberOfEdges;

        long startTime = System.currentTimeMillis();

        EvolveHelper.mainHelper(newArgs);

        long afterTime = System.currentTimeMillis();

        System.out.println("Time taken in millis: " + (afterTime - startTime));

        // #1 return result subgraph in readable form
        List<Boolean> bestIndividual = getBestIndividual(statFile);
        writeResults(bestIndividual);

        // TODO:
        // int heaviestComponentNumber = (int) graph.getHeaviestComponentInfo(bestIndividual)[1];
        // List<Boolean> heaviestComponent = graph.getComponentByNumber(bestIndividual, heaviestComponentNumber);
        // writeResults(heaviestComponent);

        // #2 write bestIndividual components to file
        List<Integer> edgesToComponents = graph.safeFindComponents(bestIndividual);
        writeConnectedComponents(edgesToComponents, bestIndividual);

        // #3 (stat by first objective)
        getNumberToFitnessStat();
    }

    private void writeConnectedComponents(List<Integer> edgesToComponents, List<Boolean> mask) {
        Map<Integer, List<Integer>> components = new HashMap<>();
        for (int i = 0; i < edgesToComponents.size(); ++i) {
            if (edgesToComponents.get(i) != 0) {
                int key = edgesToComponents.get(i);
                if (!components.containsKey(key)) {
                    components.put(key, new ArrayList<Integer>());
                }
                components.get(key).add(i);
            }
        }

        try (FileOutputStream fos   = new FileOutputStream(SOLUTION_COMPONENTS_FILE);
             BufferedWriter bw      = new BufferedWriter(new OutputStreamWriter(fos))) {

            for (Map.Entry<Integer, List<Integer>> entry : components.entrySet()) {
                int key = entry.getKey();
                double fitness = graph.getFitnessOfComponent(mask, edgesToComponents, key);
                List<Integer> edges = entry.getValue();
                bw.write("Key: " + key + ". Fitness: " + fitness + ".\n");
                for (int i = 0; i < edges.size(); ++i) {
                    Edge<Integer> edge = (Edge) graph.edges.get(edges.get(i));
                    bw.write(edge.first + " " + edge.second + "\n");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
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

    private List<Boolean> getBestIndividual(String fileName) {
        String lastLine = "";
        try (FileInputStream fs = new FileInputStream(fileName)) {
            BufferedReader br = new BufferedReader(new InputStreamReader(fs));
            String line;

            while ((line = br.readLine()) != null) {
                lastLine = line;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        if (!"".equals(lastLine)) {
            List<Boolean> mask = new ArrayList<>(lastLine.length());
            for (int i = 0; i < lastLine.length(); ++i) {
                if (lastLine.charAt(i) == '0') {
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
        Map<String, String> verticesToSignals = graph.getVerticesToSignalsWithNaN(mask);
        try (FileOutputStream fos   = new FileOutputStream(VERTICES_OUTFILE_NAME);
             BufferedWriter bw      = new BufferedWriter(new OutputStreamWriter(fos))) {

            bw.write("#label" + "\t" + "score1");
            bw.newLine();
            for (Map.Entry<String, String> entry : verticesToSignals.entrySet()) {
                bw.write(entry.getKey() + "\t" + graph.signals.get(entry.getValue()));
                bw.newLine();
            }

            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // edges
        List<Pair<Edge<String>, String>> edgesToSignals = graph.getEdgesToSignalsWithNaN(mask);
        try (FileOutputStream fos   = new FileOutputStream(EDGES_OUTFILE_NAME);
             BufferedWriter bw      = new BufferedWriter(new OutputStreamWriter(fos))) {

            bw.write("#nodeA" + "\t" + "nodeB" + "\t" + "score1");
            bw.newLine();
            for (int i = 0; i < edgesToSignals.size(); i++) {
                Pair<Edge<String>, String> edge = edgesToSignals.get(i);
                bw.write(edge.first.first + "\t" + edge.first.second + "\t" + graph.signals.get(edge.second));
                bw.newLine();
            }

            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Pattern MY_PATTERN = Pattern.compile("(\\[)(.*?)((\\])|( ))");
    public void getNumberToFitnessStat() {
        try (BufferedReader br  = new BufferedReader(new FileReader(statFile));
             PrintWriter out    = new PrintWriter(new File(NUMBER_TO_FITNESS_CSV))) {

            String line;
            int counter = 1;
            while ((line = br.readLine()) != null) {
                Matcher m = MY_PATTERN.matcher(line);
                while (m.find()) {
                    String s = m.group(2);
                    double d = Double.parseDouble(s);
                    int i = (int) d;
                    out.println(counter + "," + i);
                    ++counter;
                }
            }

        } catch (IOException e) {
            System.out.println("Error occurred while reading from file " + statFile + ".");
            e.printStackTrace();
        }

    }
}
