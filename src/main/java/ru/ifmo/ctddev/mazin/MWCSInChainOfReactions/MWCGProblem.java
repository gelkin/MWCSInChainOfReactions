package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import ec.EvolutionState;
import ec.Individual;
import ec.Problem;
import ec.multiobjective.MultiObjectiveFitness;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import ec.vector.BitVectorIndividual;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.HashMap;

public class MWCGProblem extends Problem implements SimpleProblemForm {
    private Graph graph;
    private static final String VERTICES_FILE = "./src/main/resources/nodes.txt";
    private static final String EDGES_FILE = "./src/main/resources/edges.txt";

    @Override
    public void setup(final EvolutionState state, final Parameter base) {
        // one-thread version
        initGraphFromFiles(state);
    }

    private void initGraphFromFiles(final EvolutionState state) {
        // Vertices first:
        Map<Signal, Double> signals = new HashMap<>();
        
        try {
            System.out.println((new File(".")).getCanonicalPath());
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        LinkedHashMap<String, Signal> verticesToSignal = new LinkedHashMap<>();
        try(BufferedReader br = new BufferedReader(new FileReader(VERTICES_FILE))) {
            String line = br.readLine();

            while (line != null) {
                String[] vertexInfo = line.split(" ");
                String vertex = vertexInfo[0];
                double weight = Double.parseDouble(vertexInfo[1]);

                Signal signal = new Signal(weight);
                signals.put(signal, weight);

                verticesToSignal.put(vertex, signal);

                line = br.readLine();
            }

        } catch (FileNotFoundException e) {
            state.output.fatal("Cannot find file " + VERTICES_FILE + ".", null);
            e.printStackTrace();
        } catch (IOException e) {
            state.output.fatal("Error occurred while reading from file " + VERTICES_FILE + ".", null);
            e.printStackTrace();
        }

        // Edges
        LinkedHashMap<String, LinkedHashMap<String, Signal>> edgesAsMap = new LinkedHashMap<>();
        List<Edge<String>> edges = new ArrayList<>();
        Map<String, Map<String, Integer>> edgesToIndex = new HashMap<>();

        int edgeNumber = 0;
        try(BufferedReader br = new BufferedReader(new FileReader(EDGES_FILE))) {
            String line = br.readLine();

            while (line != null) {
                String[] edgeInfo = line.split(" ");
                String from = edgeInfo[0];
                String to = edgeInfo[1];
                double weight = Double.parseDouble(edgeInfo[2]);
                Signal signal = new Signal(weight);
                signals.put(signal, weight);

                if (!edgesAsMap.containsKey(from)) {
                    edgesAsMap.put(from, new LinkedHashMap<String, Signal>());
                    edgesToIndex.put(from, new HashMap<String, Integer>());
                }

                if (!edgesAsMap.containsKey(to)) {
                    edgesAsMap.put(to, new LinkedHashMap<String, Signal>());
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
            state.output.fatal("Cannot find file " + EDGES_FILE + ".", null);
            e.printStackTrace();
        } catch (IOException e) {
            state.output.fatal("Error occurred while reading from file " + EDGES_FILE + ".", null);
            e.printStackTrace();
        }

        graph = new Graph<>(edgesAsMap, verticesToSignal, signals, edges, edgesToIndex);
    }

    @Override
    public void evaluate(final EvolutionState state,
                         final Individual individual,
                         final int subPopulation,
                         final int threadNum) {

        if (individual.evaluated) {
            return;
        }

        if (!(individual instanceof BitVectorIndividual)) {
            state.output.fatal("Whoa!  It's not a BitVectorIndividual!!!", null);
        }

        BitVectorIndividual ind2 = (BitVectorIndividual) individual;

        List<Boolean> ind2AsList = new ArrayList(ind2.genome.length);
        for (int i = 0; i < ind2.genome.length; ++i) {
            ind2AsList.add(ind2.genome[i]);
        }

        if (!(ind2.fitness instanceof MultiObjectiveFitness)) {
            state.output.fatal("Whoa!  It's not a MultiObjectiveFitness!!!", null);
        }

        double[] objectives = graph.multiObjectiveFitness(ind2AsList);

        ((MultiObjectiveFitness) ind2.fitness).setObjectives(state, objectives);

        ind2.evaluated = true;
    }

}
