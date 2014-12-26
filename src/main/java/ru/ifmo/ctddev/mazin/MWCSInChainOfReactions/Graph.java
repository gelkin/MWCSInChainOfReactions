package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import java.util.*;

public class Graph<V> implements Cloneable {
    private final LinkedHashMap<V, LinkedHashMap<V, Signal> > edgesAsMap;
    private final LinkedHashMap<V, Signal> verticesToSignals;
    private final Map<Signal, Double> signals;
    private final List<Edge<V>> edges;
    private final Map<V, Map<V, Integer>> edgesToIndex;


    private static final double SPECIAL_VERTEX_WEIGHT = -20.0;

    public Graph(LinkedHashMap<V, LinkedHashMap<V, Signal>> edgesAsMap,
                 LinkedHashMap<V, Signal> verticesToSignals,
                 Map<Signal, Double> signals,
                 List<Edge<V>> edges,
                 Map<V, Map<V, Integer>> edgesToIndex) {

        this.edgesAsMap = edgesAsMap;
        this.verticesToSignals = verticesToSignals;
        this.signals = signals;
        this.edges = edges;
        this.edgesToIndex = edgesToIndex;
    }

    private List<Boolean> mask;
    private List<Integer> edgeToComponent;
    private int compNumber;
    private int majorComponentIndex = -1;
    // returns weight of connected component max-weight
    public double fitness(List<Boolean> mask) {
        this.mask = mask;

        findComponents();

        int[] freq = new int[compNumber + 1];
        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                freq[edgeToComponent.get(i)]++;
            }
        }
        
        int biggestCompNumber = -1;
        int edgeNumber = 0;
        for (int i = 0; i < compNumber; ++i) {
            if (freq[i] > edgeNumber) {
                biggestCompNumber = i;
            }
        }
        
        double fitness = getFitnessOfComponent(biggestCompNumber);
        
        /*
        for (int i = 0; i < mask.size(); ++i) {
            if (edgeToComponent.get(i) == biggestCompNumber) {
                fitn
            }
            double newFitness = getFitnessOfComponent(i);
            if (fitness < newFitness) {
                fitness = newFitness;
                majorComponentIndex = i;
            }
        }
        */

        return fitness;
    }

    // returns array with values [fitness(), numberOfComponents, edgesNumber]
    public double[] multiObjectiveFitness(List<Boolean> mask) {
        double[] objectives = new double[3];
        objectives[2] = fitness(mask) / 100;
       
        // returns: [numberOfEdges, numberOfVertices]
    
        // objectives[1] = compNumber;
        
        int[] componentSize = new int[2];
        Map<V, Boolean> isCountedVertex = new HashMap<>();
        for (int i = 0; i < mask.size(); ++i) {
            if (edgeToComponent.get(i) == majorComponentIndex) {
                ++componentSize[0];
                V first = edges.get(i).first;
                V second = edges.get(i).second;

                if (!isCountedVertex.containsKey(first)) {
                    isCountedVertex.put(first, true);
                    ++componentSize[1];
                }

                if (!isCountedVertex.containsKey(second)) {
                    isCountedVertex.put(second, true);
                    ++componentSize[1];
                }
            }
        }

        // put numberOfEdges as objective
        objectives[0] = componentSize[0];
        objectives[1] = componentSize[1];

        return objectives;
    }


    // returns: [numberOfEdges, numberOfVertices]
    public int[] getMajorComponentSize(List<Boolean> mask) {
        fitness(mask);

        int[] componentSize = new int[2];
        Map<V, Boolean> isCountedVertex = new HashMap<>();
        for (int i = 0; i < mask.size(); ++i) {
            if (edgeToComponent.get(i) == majorComponentIndex) {
                ++componentSize[0];
                V first = edges.get(i).first;
                V second = edges.get(i).second;

                if (!isCountedVertex.containsKey(first)) {
                    isCountedVertex.put(first, true);
                    ++componentSize[1];
                }

                if (!isCountedVertex.containsKey(second)) {
                    isCountedVertex.put(second, true);
                    ++componentSize[1];
                }
            }
        }

        return componentSize;
    }  
    

    public double getEdgesNumber(List<Boolean> mask) {
        double edgesNumber = 0.0;
        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                edgesNumber += 1.0;
            }
        }
        return edgesNumber;
    }

    private void findComponents() {
        edgeToComponent = new ArrayList<>(Collections.nCopies(edges.size(), 0));

        compNumber = 0;
        for (int i = 0; i < edges.size(); ++i) {
            if (mask.get(i) && edgeToComponent.get(i) == 0) {
                ++compNumber;
                dfs(i, compNumber, edgeToComponent);
            }
        }
    }

    private void dfs(int edgeNumber, int compNumber, List<Integer> edgeToComponent) {
        edgeToComponent.set(edgeNumber, compNumber);

        Edge<V> edge = edges.get(edgeNumber);

        V first = edge.first;
        for (Map.Entry<V, Signal> entry : edgesAsMap.get(first).entrySet()) {
            V firstTo = entry.getKey();
            int indexOfEdge = edgesToIndex.get(first).get(firstTo);
            if (mask.get(indexOfEdge) && edgeToComponent.get(indexOfEdge) != compNumber) {
                dfs(indexOfEdge, compNumber, edgeToComponent);
            }
        }

        V second = edge.first;
        for (Map.Entry<V, Signal> entry : edgesAsMap.get(second).entrySet()) {
            V secondTo = entry.getKey();
            int indexOfEdge = edgesToIndex.get(second).get(secondTo);
            // if edge in graph and isn't counted then
            if (mask.get(indexOfEdge) && edgeToComponent.get(indexOfEdge) != compNumber) {
                dfs(indexOfEdge, compNumber, edgeToComponent);
            }
        }
    }

    private double getFitnessOfComponent(int compIndex) {
        double fitness = 0.0;
        Set<Signal> uniqueSignals = new HashSet<>();
        for (int i = 0; i < edges.size(); ++i) {
            if (mask.get(i) && edgeToComponent.get(i) == compIndex) {
                Edge<V> edge = edges.get(i);
                uniqueSignals.add(edgesAsMap.get(edge.first).get(edge.second));

                if (signals.get(verticesToSignals.get(edge.first)) != SPECIAL_VERTEX_WEIGHT) {
                    uniqueSignals.add(verticesToSignals.get(edge.first));
                } else {
                    fitness += SPECIAL_VERTEX_WEIGHT;
                }

                if (signals.get(verticesToSignals.get(edge.second)) != SPECIAL_VERTEX_WEIGHT) {
                    uniqueSignals.add(verticesToSignals.get(edge.second));
                } else {
                    fitness += SPECIAL_VERTEX_WEIGHT;
                }
            }
        }

        Iterator<Signal> it = uniqueSignals.iterator();
        while (it.hasNext()) {
            fitness += signals.get(it.next());
        }

        return fitness;
    }
}
