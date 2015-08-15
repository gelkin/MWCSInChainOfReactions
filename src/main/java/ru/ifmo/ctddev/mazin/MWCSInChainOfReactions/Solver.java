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
    private static final String INITIAL_POP_PARAM = "pop.subpop.0.file=";
    private static final String SUBPOPS_SIZE_PARAM = "pop.subpop.0.size=";
    private static final String MUTATION_PROB_PARAM = "pop.subpop.0.species.mutation-prob=";
    private static final int SUBPOPULATION_SIZE = 50;
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
        } else if (args.length > 4) {
            System.out.println("Too many arguments");
            return;
        }

        try {
            String verticesFile = args[0];
            String edgesFile = args[1];
            String signalsFile = args[2];

            graph = initGraphFromFiles(verticesFile, edgesFile, signalsFile);
            setGraph(graph);
            int numberOfEdges = getNumberOfLinesInFile(edgesFile);

            List<String> newArgs = new ArrayList<>();
            // add parameters file
            newArgs.add(FILE_PARAM);
            newArgs.add(parametersFile);
            // set number of edges in graph - individual.genome size
            newArgs.add(PARAM_IDENTIFIER);
            newArgs.add(GENOME_SIZE_PARAM + numberOfEdges);
            // set size of (sub)population
            newArgs.add(PARAM_IDENTIFIER);
            newArgs.add(SUBPOPS_SIZE_PARAM + SUBPOPULATION_SIZE);
            // set mutation probability (depends on number of edges)
            newArgs.add(PARAM_IDENTIFIER);
            newArgs.add(MUTATION_PROB_PARAM + (1.0 / numberOfEdges));
            if (args.length == 4) {
                // We have to set initial population
                String initialPopFile = args[3];
                createFileWithInitialPopInNSGA2Style(initialPopFile);

                // set initial population
                newArgs.add(PARAM_IDENTIFIER);
                newArgs.add(INITIAL_POP_PARAM + "$" + INITIAL_POP_ECJ_STYLE_FILE);
            }

            EvolveHelper.mainHelper(newArgs.toArray(new String[newArgs.size()]));

            // #1 Get best connected individual (with 1.0 connected component)
            List<Boolean> bestIndividual = getBestConnectedIndividual();
            writeResults(bestIndividual);
            // #2 Check: if we got connected individual
            List<Integer> edgesToComponents = (List) graph.findComponents(bestIndividual).first;
            writeConnectedComponents(edgesToComponents, bestIndividual);
            // #3 Write statistics by first objective
            writeNumberToFitnessStat();

        } catch (MWCSException me) {
            System.out.println("Exception caught:\n" + me.toString());
        }
    }

    private Graph initGraphFromFiles(String verticesFile, String edgesFile, String signalsFile) throws MWCSException {
        Map<String, String> verticesToSignal = initVerticesFromFile(verticesFile);
        Map<String, Double> signals = initSignalsFromFile(signalsFile);

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
            throw new MWCSException("Cannot find file " + edgesFile + ".");
        } catch (IOException e) {
            throw new MWCSException("Error occurred while reading from file " + edgesFile + ".");
        }

        return new Graph<>(edgesAsMap, verticesToSignal, signals, edges, edgesToIndex);
    }

    private Map<String, String> initVerticesFromFile(String verticesFile) throws MWCSException {
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
            throw new MWCSException("Cannot find file " + verticesFile + ".");
        } catch (IOException e) {
            throw new MWCSException("Error occurred while reading from file " + verticesFile + ".");
        }
        return verticesToSignal;
    }

    private Map<String, Double> initSignalsFromFile(String signalsFile) throws MWCSException {
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
            throw new MWCSException("Cannot find file " + signalsFile + ".");
        } catch (IOException e) {
            throw new MWCSException("Error occurred while reading from file " + signalsFile + ".");
        }
        return signals;
    }

    private void setGraph(Graph graph) {
        MWCGProblem.setGraph(graph);
        GraphMutatorPipeline.setGraph(graph);
        GraphIndividual.setGraph(graph);
    }

    private int getNumberOfLinesInFile(String fileName) throws MWCSException {
        int linesNumber = 0;
        try(BufferedReader br = new BufferedReader(new FileReader(fileName))) {
            String line = br.readLine();

            while (line != null) {
                ++linesNumber;
                line = br.readLine();
            }
        } catch (FileNotFoundException e) {
            throw new MWCSException("Cannot find file " + fileName + ".");
        } catch (IOException e) {
            throw new MWCSException("Error occurred while reading from file.");
        }

        return linesNumber;
    }

    private void createFileWithInitialPopInNSGA2Style(String initialPopFile) throws MWCSException {
        Map<String, Map<String, Boolean>> initialPopAsEdges = initInitialPopAsEdges(initialPopFile);

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

            bw.write("Number of Individuals: i" + SUBPOPULATION_SIZE + "|");
            bw.newLine();
            for (int i = 0; i < SUBPOPULATION_SIZE; ++i) {
                bw.write("Individual Number: i" + i + "|");
                bw.newLine();
                bw.write("Evaluated: T");
                bw.newLine();
                bw.write("Fitness: [d|" + fitness[0] + "|d|" + fitness[1] + "|]");
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
            throw new MWCSException("Error occurred while writing edges to file " + INITIAL_POP_ECJ_STYLE_FILE + ".");
        }
    }

    private Map<String, Map<String, Boolean>> initInitialPopAsEdges(String initialPopFile) throws MWCSException {
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
            throw new MWCSException("Cannot find file " + initialPopFile + ".");
        } catch (IOException e) {
            throw new MWCSException("Error occurred while reading from file " + initialPopFile + ".");
        }
        return edgesAsMap;
    }

    private void writeResults(List<Boolean> mask) throws MWCSException {
        writeResultVertices(mask);
        writeResultEdges(mask);
    }

    private void writeResultVertices(List<Boolean> mask) throws MWCSException {
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
            throw new MWCSException("Error occurred while writing vertices to file " + VERTICES_OUTFILE_NAME + ".");
        }
    }

    private void writeResultEdges(List<Boolean> mask) throws MWCSException {
        List<Pair<Edge<String>, String>> edgesToSignals = graph.getEdgesToSignalsWithNaN(mask);
        try (FileOutputStream fos   = new FileOutputStream(EDGES_OUTFILE_NAME);
             BufferedWriter bw      = new BufferedWriter(new OutputStreamWriter(fos))) {

            bw.write("#nodeA" + "\t" + "nodeB" + "\t" + "score1");
            bw.newLine();
            for (Pair<Edge<String>, String> edge : edgesToSignals) {
                bw.write(edge.first.first + "\t" + edge.first.second + "\t" + graph.signals.get(edge.second));
                bw.newLine();
            }

            bw.close();
        } catch (IOException e) {
            throw new MWCSException("Error occurred while writing edges to file " + EDGES_OUTFILE_NAME + ".");
        }
    }

    private static Pattern ALL_STAT_PATTERN = Pattern.compile("\\[(.*?) (.*?)\\]");
    private List<Boolean> getBestConnectedIndividual() throws MWCSException {
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
            throw new MWCSException("Error occurred while reading from file " + statFile + ".");
        }
    }

    private void writeConnectedComponents(List<Integer> edgesToComponents, List<Boolean> mask) throws MWCSException {
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
            throw new MWCSException("Error occurred while writing to file " + RESULT_COMPONENTS_FILE + ".");
        }
    }

    private static Pattern NUMBER_TO_FITNESS_PATTERN = Pattern.compile("(\\[)(.*?)((\\])|( ))");
    public void writeNumberToFitnessStat() throws MWCSException {
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
            throw new MWCSException("Error occurred while reading from file " + statFile + ".");
        }
    }
}
