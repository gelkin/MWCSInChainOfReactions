package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import java.util.*;

/**
 *
 * @param <V> - vertex type
 * @param <S> - signal type
 */
public class Graph<V, S> implements Cloneable {
    private final LinkedHashMap<V, LinkedHashMap<V, S> > edgesAsMap;
    private final LinkedHashMap<V, S> verticesToSignals;
    private final Map<S, Double> signals;
    private final List<Edge<V>> edges;
    private final Map<V, Map<V, Integer>> edgesToIndex;
    private final List<Integer> edgeToComponentInWholeGraph;
    private final int componentsNumberInWholeGraph;
    private int[] orderToNumberForComponents;

    private static final double SPECIAL_VERTEX_WEIGHT = -20.0;

    public Graph(LinkedHashMap<V, LinkedHashMap<V, S>> edgesAsMap,
                 LinkedHashMap<V, S> verticesToSignals,
                 Map<S, Double> signals,
                 List<Edge<V>> edges,
                 Map<V, Map<V, Integer>> edgesToIndex) {

        this.edgesAsMap = edgesAsMap;
        this.verticesToSignals = verticesToSignals;
        this.signals = signals;
        this.edges = edges;
        this.edgesToIndex = edgesToIndex;

        edgeToComponentInWholeGraph = findComponents();
        componentsNumberInWholeGraph = componentsNumber;
        orderToNumberForComponents = orderToNumber();
    }

    /*
     * Main fitness function.
     * 
     * @param 'p':
     * p < 0 -  return  DoubleMinValue
     * p == 0 - search in whole graph
     * 0 < p < n - search in first 'p' connected components
     * p >= n - search in whole graph
    */
    public double fitness(List<Boolean> originalMask, int p) {
        if (p < 0) {
            return Double.MIN_VALUE;
        }

        List<Boolean> mask = updateMask(originalMask, p);
        List<Integer> edgeToComponent = findComponents(mask);
        int biggestComponentNumber = getBiggestByEdgeComponentNumber(mask, edgeToComponent);
        double biggestComponentFitness = getFitnessOfComponent(biggestComponentNumber, mask, edgeToComponent);

        return biggestComponentFitness;
    }

    public List<Boolean> getBiggestComponent(List<Boolean> mask) {
        List<Integer> edgeToComponent = findComponents(mask);
        int biggestcomponentNumber = getBiggestByEdgeComponentNumber(mask, edgeToComponent);

        List<Boolean> biggestComponentMask = new ArrayList<>(mask.size());
        for (int i = 0; i < edges.size(); ++i) {
            if (edgeToComponent.get(i) == biggestcomponentNumber) {
                biggestComponentMask.add(true);
            } else {
                biggestComponentMask.add(false);
            }
        }

        return biggestComponentMask;
    }

    // todo mb remove
    // returns array with values [fitness(), numberOfComponents, edgesNumber]
    private double[] multiObjectiveFitness(List<Boolean> mask,
                                          List<Integer> edgeToComponent) {

        double[] objectives = new double[3];
        objectives[2] = fitness(mask, 0) / 100;

        // returns: [numberOfEdges, numberOfVertices]
    
        // objectives[1] = componentNumber;
        
        int[] componentSize = new int[2];
        Map<V, Boolean> isCountedVertex = new HashMap<>();
        for (int i = 0; i < mask.size(); ++i) {
            if (edgeToComponent.get(i) == /* majorComponentIndex todo*/0) {
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

    private int getBiggestByEdgeComponentNumber(List<Boolean> mask,
                                                List<Integer> edgeToComponent) {
        int[] compToEdgesNumber = new int[edges.size()];
        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                ++compToEdgesNumber[edgeToComponent.get(i)];
            }
        }

        int res = -1;
        int edgesNumber = 0;
        for (int i = 0; i < edges.size(); ++i) {
            if (compToEdgesNumber[i] > edgesNumber) {
                res = i;
                edgesNumber = compToEdgesNumber[i];
            }
        }

        return res;
    }

    private double getHeaviestComponentFitness(List<Boolean> mask,
                                               List<Integer> edgeToComponent) {
        double fitness = Double.MIN_VALUE;
        for (int i = 0; i < componentsNumber; ++i) {
            double newFitness = getFitnessOfComponent(i, mask, edgeToComponent);
            if (fitness < newFitness) {
                fitness = newFitness;
            }
        }

        return fitness;

    }

    // returns: [numberOfEdges, numberOfVertices]
    private int[] getBiggestByEdgeComponentSize(List<Boolean> mask,
                                         List<Integer> edgeToComponent) {
        fitness(mask, 0);

        int[] componentSize = new int[2];
        Map<V, Boolean> isCountedVertex = new HashMap<>();
        for (int i = 0; i < mask.size(); ++i) {
            if (edgeToComponent.get(i) == /*majorComponentIndex*/ 0) {
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
    

    public int getEdgesNumber() {
        return edges.size();
    }

    private List<Integer> findComponents() {
        List<Boolean> wholeGraphMask = new ArrayList<>(Collections.nCopies(edges.size(), true));
        return findComponents(wholeGraphMask);
    }

    private int componentsNumber;
    // return: list.get(i) = number_of_component , if edge #i match the 'mask'
    //                       0 , if edge #i doesn't match the 'mask'
    private List<Integer> findComponents(List<Boolean> mask) {
        List<Integer> edgeToComponent = new ArrayList<>(Collections.nCopies(edges.size(), 0));

        int componentNumber = 0;
        for (int i = 0; i < edges.size(); ++i) {
            if (mask.get(i) && edgeToComponent.get(i) == 0) {
                ++componentNumber;
                dfs(i, componentNumber, mask, edgeToComponent);
            }
        }
        
        componentsNumber = componentNumber;

        return edgeToComponent;
    }

    private void dfs(int edgeNumber, int componentNumber, List<Boolean> mask, List<Integer> edgeToComponent) {
        edgeToComponent.set(edgeNumber, componentNumber);

        Edge<V> edge = edges.get(edgeNumber);

        V first = edge.first;
        for (Map.Entry<V, S> entry : edgesAsMap.get(first).entrySet()) {
            V firstTo = entry.getKey();
            int indexOfEdge = edgesToIndex.get(first).get(firstTo);
            if (mask.get(indexOfEdge) && edgeToComponent.get(indexOfEdge) != componentNumber) {
                dfs(indexOfEdge, componentNumber, mask, edgeToComponent);
            }
        }

        V second = edge.first;
        for (Map.Entry<V, S> entry : edgesAsMap.get(second).entrySet()) {
            V secondTo = entry.getKey();
            int indexOfEdge = edgesToIndex.get(second).get(secondTo);
            // if edge in graph and isn't counted then
            if (mask.get(indexOfEdge) && edgeToComponent.get(indexOfEdge) != componentNumber) {
                dfs(indexOfEdge, componentNumber, mask, edgeToComponent);
            }
        }
    }

    private double getFitnessOfComponent(int compIndex, List<Boolean> mask, List<Integer> edgeToComponent) {
        double fitness = 0.0;
        Set<S> uniqueSignals = new HashSet<>();
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

        Iterator<S> it = uniqueSignals.iterator();
        while (it.hasNext()) {
            fitness += signals.get(it.next());
        }

        return fitness;
    }

    // Returns mask with edges only in first 'p' largest connected components.
    private List<Boolean> updateMask(List<Boolean> mask, int p) {
        if (p == 0) {
            return mask;
        }

        List<Boolean> newMask = new ArrayList<>(Collections.nCopies(mask.size(), false));
        mainLoop:
        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                // todo optimize number of operations
                for (int j = 0; j < p; ++j) {
                    if (orderToNumberForComponents[j] == (edgeToComponentInWholeGraph.get(i) - 1)) {
                        newMask.set(i, true);
                        continue mainLoop;
                    }
                }
            }
        }

        return newMask;
    }

    // todo down
    // Get info about subgraph presented by mask:
    public List<Pair<V, S>> getVerticesToSignalsByMask(List<Boolean> mask) {
        Set<V> vertices = new HashSet<>();
        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                vertices.add(edges.get(i).first);
                vertices.add(edges.get(i).second);
            }
        }

        List<Pair<V, S>> result = new ArrayList<>(vertices.size());
        Iterator<V> it = vertices.iterator();
        while (it.hasNext()) {
            V vertex = it.next();
            result.add(new Pair<V, S>(vertex, verticesToSignals.get(vertex)));
        }

        return result;
    }

    public List<Pair<Edge<V>, S>> getEdgesToSignalsByMask(List<Boolean> mask) {
        List<Pair<Edge<V>, S>> result = new ArrayList<>();

        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                Edge e = edges.get(i);
                result.add(new Pair(e, edgesAsMap.get(e.first).get(e.second)));
            }
        }

        return result;
    }

    public List<Pair<S, Double>> getSignalsToValuesByMask(List<Boolean> mask) {
        Set<S> maskSignals = new HashSet<>();
        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                maskSignals.add(verticesToSignals.get(edges.get(i).first));
                maskSignals.add(verticesToSignals.get(edges.get(i).second));
                maskSignals.add(edgesAsMap.get(edges.get(i).first).get(edges.get(i).second));
            }
        }

        List<Pair<S, Double>> result = new ArrayList<>(maskSignals.size());
        Iterator<S> it = maskSignals.iterator();
        while (it.hasNext()) {
            S signal = it.next();
            result.add(new Pair<>(signal, signals.get(signal)));
        }

        return result;
    }
    // todo up

    // Returns array of the form:
    // <order_of_connected_component_by_size> -> <number_of_connected_component>
    private int[] orderToNumber() {
        // componentNumberToSize:
        List<ComponentInfo> components = new ArrayList<>(componentsNumberInWholeGraph);
        for (int i = 0; i < componentsNumberInWholeGraph; ++i) {
            components.add(i, new ComponentInfo(i, 0));
        }

        for (int i = 0; i < edges.size(); ++i) {
            int componentNumber = edgeToComponentInWholeGraph.get(edgesToIndex.get(edges.get(i).first).get(edges.get(i).second));
            components.get(componentNumber - 1).size++;
        }

        Collections.sort(components);

        int[] orderToNumber = new int[componentsNumberInWholeGraph];
        for (int i = 0; i < componentsNumberInWholeGraph; ++i) {
            orderToNumber[i] = components.get(i).number;
        }

        return orderToNumber;
    }

    private class ComponentInfo implements Comparable<ComponentInfo> {
        public int number;
        public int size;

        private ComponentInfo(int number, int size) {
            this.number = number;
            this.size = size;
        }

        public int compareTo(ComponentInfo element) {
            if (size > element.size) {
                return -1;
            } else if (size < element.size) {
                return 1;
            }
            return 0;
        }
    }

}
