package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import java.io.*;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Solver {
    // TODO src/main/resources/true-solution-edges.txt
    private static final String FILE_PARAM = "-file";
    private static final String PARAM_IDENTIFIER = "-p";
    private static final String GENOME_SIZE_PARAM = "pop.subpop.0.species.genome-size=";
    private static final String INITIAL_POP_PARAM = "pop.subpop.0.file=";
    private static final String SUPPOPS_SIZE_PARAM = "pop.subpop.0.size =";
    private static int SUPPOPS_SIZE = 50;
    private static final String VERTICES_OUTFILE_NAME = "./src/main/resources/result-nodes.txt";
    private static final String EDGES_OUTFILE_NAME = "./src/main/resources/result-edges.txt";
    private static final String NUMBER_TO_FITNESS_CSV = "./src/main/resources/number_to_fitness.csv";
    private static final String RESULT_COMPONENTS_FILE = "./src/main/resources/result-components.stat";
    private static final String INITIAL_POP_ECJ_STYLE_FILE = "./src/main/resources/init-pop-ecj-style.txt";

    private static String parametersFile = "./src/main/resources/mwcs-multi.params";
    private static String statFile = "./src/main/resources/mwcs-multi.stat";

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
        if (graph == null) {
            return;
        }

        setGraph(graph);

        int numberOfEdges = getNumberOfLinesInFile(edgesFile);
        if (numberOfEdges == -1) {
            return;
        }

        final String[] newArgs;
        if (args.length == 3) {
            newArgs = new String[6];
        } else {
            // We have a file to set initial population
            newArgs = new String[8];

            String initialPopFile = args[3];
            if (!createFileWithInitialPop(initialPopFile)) {
                return;
            }
            newArgs[6] = PARAM_IDENTIFIER;
            newArgs[7] = INITIAL_POP_PARAM + "$" + INITIAL_POP_ECJ_STYLE_FILE;
        }
        newArgs[0] = FILE_PARAM;
        newArgs[1] = parametersFile;
        newArgs[2] = PARAM_IDENTIFIER;
        newArgs[3] = GENOME_SIZE_PARAM + numberOfEdges;
        newArgs[4] = PARAM_IDENTIFIER;
        newArgs[5] = SUPPOPS_SIZE_PARAM + SUPPOPS_SIZE;

        EvolveHelper.mainHelper(newArgs);

        // #1 Get best connected individual (with 1.0 connected component)
        List<Boolean> bestIndividual = getBestConnectedIndividual();
        if (bestIndividual == null || !writeResults(bestIndividual)) {
            return;
        }
        // #2 Check: if we got connected individual
        List<Integer> edgesToComponents = (List) graph.findComponents(bestIndividual).first;
        if (!writeConnectedComponents(edgesToComponents, bestIndividual)) {
            return;
        }
        // #3 Write statistics by first objective
        writeNumberToFitnessStat();
    }

    private Graph initGraphFromFiles(String verticesFile, String edgesFile, String signalsFile) {
        Map<String, String> verticesToSignal = initVerticesFromFile(verticesFile);
        if (verticesToSignal == null) {
            return null;
        }
        Map<String, Double> signals = initSignalsFromFile(signalsFile);
        if (signals == null) {
            return null;
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
            return null;
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file " + edgesFile + ".");
            e.printStackTrace();
            return null;
        }

        return new Graph<>(edgesAsMap, verticesToSignal, signals, edges, edgesToIndex);
    }

    private Map<String, String> initVerticesFromFile(String verticesFile) {
        Map<String, String> verticesToSignal = new HashMap<>();
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
            return null;
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file " + verticesFile + ".");
            e.printStackTrace();
            return null;
        }
        return verticesToSignal;
    }

    private Map<String, Double> initSignalsFromFile(String signalsFile) {
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
            System.out.println("Cannot find file " + signalsFile + ".");
            e.printStackTrace();
            return null;
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file " + signalsFile + ".");
            e.printStackTrace();
            return null;
        }
        return signals;
    }

    private boolean createFileWithInitialPop(String initialPopFile) {
        Map<String, Map<String, Boolean>> initialPopAsEdges = initInitialPopAsEdges(initialPopFile);
        if (initialPopAsEdges == null) {
            return false;
        }

        List<Boolean> mask = new ArrayList<>(graph.edges.size());
        StringBuilder ecjStyleIndividual = new StringBuilder(graph.edges.size());
        ecjStyleIndividual.append("i" + graph.edges.size() + "|");
        for (Edge<String> edge : (List<Edge<String>>) graph.edges) {
            if (initialPopAsEdges.get(edge.first).get(edge.second)) {
                mask.add(true);
                ecjStyleIndividual.append("i1|");
            } else {
                mask.add(false);
                ecjStyleIndividual.append("i0|");
            }
        }

        double[] fitness = graph.fitness(mask, 0);

        // Need initial population to have the same ordering as in whole graph
        try (FileOutputStream fos   = new FileOutputStream(INITIAL_POP_ECJ_STYLE_FILE);
             BufferedWriter bw      = new BufferedWriter(new OutputStreamWriter(fos))) {

            bw.write("Number of Individuals: i" + SUPPOPS_SIZE + "|");
            bw.newLine();
            for (int i = 0; i < SUPPOPS_SIZE; ++i) {
                bw.write("Individual Number: i" + i + "|");
                bw.newLine();
                bw.write("Evaluated: T");
                bw.newLine();
                bw.write("Fitness: [d|" + fitness[0] + "| d|" + fitness[1] + "|]");
                bw.newLine();
                bw.write("Rank: i0|");
                bw.newLine();
                bw.write("Sparsity:  d|Infinity|");
                bw.newLine();
                bw.write(ecjStyleIndividual.toString());
                bw.newLine();
            }

            bw.close();
        } catch (IOException e) {
            System.out.println("Error occurred while writing edges to file " + INITIAL_POP_ECJ_STYLE_FILE + ".");
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private Map<String, Map<String, Boolean>> initInitialPopAsEdges(String initialPopFile) {
        Map<String, Map<String, Boolean>> edgesAsMap = new HashMap<>();
        try(BufferedReader br = new BufferedReader(new FileReader(initialPopFile))) {
            String line = br.readLine();

            while (line != null) {
                String[] edgeInfo = line.split("\\t");
                String from = edgeInfo[0];
                String to = edgeInfo[1];
                String signal = edgeInfo[2];

                if (!edgesAsMap.containsKey(from)) {
                    edgesAsMap.put(from, new HashMap<String, Boolean>());
                }

                if (!edgesAsMap.containsKey(to)) {
                    edgesAsMap.put(to, new HashMap<String, Boolean>());
                }

                if (signal.equals("NaN")) {
                    edgesAsMap.get(from).put(to, false);
                    edgesAsMap.get(to).put(from, false);
                } else {
                    edgesAsMap.get(from).put(to, true);
                    edgesAsMap.get(to).put(from, true);
                }

                line = br.readLine();
            }
        } catch (FileNotFoundException e) {
            System.out.println("Cannot find file " + initialPopFile + ".");
            e.printStackTrace();
            return null;
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file " + initialPopFile + ".");
            e.printStackTrace();
            return null;
        }
        return edgesAsMap;
    }

    private void setGraph(Graph graph) {
        MWCGProblem.setGraph(graph);
        GraphMutatorPipeline.setGraph(graph);
        GraphIndividual.setGraph(graph);
    }

    private int getNumberOfLinesInFile(String fileName) {
        int linesNumber = 0;
        try(BufferedReader br = new BufferedReader(new FileReader(fileName))) {
            String line = br.readLine();

            while (line != null) {
                ++linesNumber;
                line = br.readLine();
            }
        } catch (FileNotFoundException e) {
            System.out.println("Cannot find file " + fileName + ".");
            e.printStackTrace();
            linesNumber = -1;
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file.");
            e.printStackTrace();
            linesNumber = -1;
        }

        return linesNumber;
    }

    private boolean writeResults(List<Boolean> mask) {
        if (!writeResultVertices(mask)) {
            return false;
        }
        return writeResultEdges(mask);
    }

    private boolean writeResultVertices(List<Boolean> mask) {
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
            System.out.println("Error occurred while writing vertices to file " + VERTICES_OUTFILE_NAME + ".");
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private boolean writeResultEdges(List<Boolean> mask) {
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
            System.out.println("Error occurred while writing edges to file " + EDGES_OUTFILE_NAME + ".");
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static Pattern ALL_STAT_PATTERN = Pattern.compile("\\[(.*?) (.*?)\\]");
    private List<Boolean> getBestConnectedIndividual() {
        try (BufferedReader br  = new BufferedReader(new FileReader(statFile))) {
            String line;
            Pair<Double, String> maxWeightSubgraph = new Pair<>(-Double.MAX_VALUE, "");
            while ((line = br.readLine()) != null) {
                String[] parts = line.split(" ");
                if ("Fitness:".equals(parts[0])) {
                    Matcher m = ALL_STAT_PATTERN.matcher(line);
                    while (m.find()) {
                        String s = m.group(1);
                        double numberOfComponents = Double.parseDouble(s);
                        s = m.group(2);
                        double weight = Double.parseDouble(s);
                        if (numberOfComponents == 1.0 && weight > maxWeightSubgraph.first) {
                            maxWeightSubgraph.first = weight;
                            br.readLine();
                            br.readLine();
                            maxWeightSubgraph.second = br.readLine();
                        }
                    }
                }
            }

            List<Boolean> mask = new ArrayList<>(maxWeightSubgraph.second.length());
            for (int i = 0; i < maxWeightSubgraph.second.length(); ++i) {
                if ('0' == maxWeightSubgraph.second.charAt(i)) {
                    mask.add(false);
                } else {
                    mask.add(true);
                }
            }

            writeConnectedComponents((List) graph.findComponents(mask).first, mask);

            return mask;
        } catch (IOException e) {
            System.out.println("Error occurred while reading from file " + statFile + ".");
            e.printStackTrace();
            return null;
        }
    }

    private boolean writeConnectedComponents(List<Integer> edgesToComponents, List<Boolean> mask) {
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

        try (FileOutputStream fos   = new FileOutputStream(RESULT_COMPONENTS_FILE);
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
            System.out.println("Error occurred while writing to file " + RESULT_COMPONENTS_FILE + ".");
            e.printStackTrace();
            return false;
        }
        return true;
    }

    private static Pattern NUMBER_TO_FITNESS_PATTERN = Pattern.compile("(\\[)(.*?)((\\])|( ))");
    public boolean writeNumberToFitnessStat() {
        try (BufferedReader br  = new BufferedReader(new FileReader(statFile));
             PrintWriter out    = new PrintWriter(new File(NUMBER_TO_FITNESS_CSV))) {

            String line;
            int counter = 1;
            while ((line = br.readLine()) != null) {
                Matcher m = NUMBER_TO_FITNESS_PATTERN.matcher(line);
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
            return false;
        }
        return true;
    }
}
