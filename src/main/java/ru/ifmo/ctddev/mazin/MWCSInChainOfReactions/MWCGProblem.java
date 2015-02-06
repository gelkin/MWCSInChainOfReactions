package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import ec.EvolutionState;
import ec.Individual;
import ec.Problem;
import ec.multiobjective.MultiObjectiveFitness;
import ec.simple.SimpleFitness;
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
    private static String verticesFile;
    private static String edgesFile;
    private static String signalsFile;

    @Override
    public void setup(final EvolutionState state, final Parameter base) {
        // one-thread version
        initGraphFromFiles(state);
    }

    private static boolean areSourceFileSet = false;
    public static void setSourceFiles(String verticesFile, String edgesFile, String signalsFile) {
        if (!areSourceFileSet) {
            MWCGProblem.verticesFile = verticesFile;
            MWCGProblem.edgesFile = edgesFile;
            MWCGProblem.signalsFile = signalsFile;
            areSourceFileSet = true;
        }
    }

    private void initGraphFromFiles(final EvolutionState state) {

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
            state.output.fatal("Cannot find file " + verticesFile + ".", null);
            e.printStackTrace();
        } catch (IOException e) {
            state.output.fatal("Error occurred while reading from file " + verticesFile + ".", null);
            e.printStackTrace();
        }

        // Edges:
        LinkedHashMap<String, LinkedHashMap<String, String>> edgesAsMap = new LinkedHashMap<>();
        List<Edge<String>> edges = new ArrayList<>();
        Map<String, Map<String, Integer>> edgesToIndex = new HashMap<>();

        System.out.println("edgesFile: " + edgesFile);

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
            state.output.fatal("Cannot find file " + edgesFile + ".", null);
            e.printStackTrace();
        } catch (IOException e) {
            state.output.fatal("Error occurred while reading from file " + edgesFile + ".", null);
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
            state.output.fatal("Cannot find file " + verticesFile + ".", null);
            e.printStackTrace();
        } catch (IOException e) {
            state.output.fatal("Error occurred while reading from file " + verticesFile + ".", null);
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

        /*
        if (!(ind2.fitness instanceof MultiObjectiveFitness)) {
            state.output.fatal("Whoa!  It's not a MultiObjectiveFitness!!!", null);
        }

        double[] objectives = graph.multiObjectiveFitness(ind2AsList);

        ((MultiObjectiveFitness) ind2.fitness).setObjectives(state, objectives);
        */

        if (!(ind2.fitness instanceof SimpleFitness)) {
            state.output.fatal("Whoa!  It's not a SimpleFitness!!!", null);
        }

        ((SimpleFitness)ind2.fitness).setFitness(state,
                                                 graph.fitness(ind2AsList, 4),
                                                 false);

        ind2.evaluated = true;
    }

}
