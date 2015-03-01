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
     * p < 0 -  return  DoubleNegativeInfinity
     * p == 0 - search in whole graph
     * 0 < p < n - search in first 'p' connected components
     * p >= n - search in whole graph
     *  where 'n' - number of connected components in graph
    */
    public double[] fitness(List<Boolean> originalMask, int p) {
        if (p < 0) {
            // Multi
            double[] result = {Double.NEGATIVE_INFINITY, 0, 0};
            return result;
        }

        List<Boolean> mask = updateMask(originalMask, p);
        List<Integer> edgeToComponent = findComponents(mask);

        // Multi:
        return multiObjectiveFitness(mask, edgeToComponent);
    }

    // returns array with values [fitness(), edgesOfComponent, verticesOfComponent]
    private static final int NUM_OBJECTIVES = 1;
    private double[] multiObjectiveFitness(List<Boolean> mask,
                                           List<Integer> edgeToComponent) {
        double[] objectives = new double[NUM_OBJECTIVES];

        // double[] info = getHeaviestComponentInfo(mask, edgeToComponent);

        int x = getBiggestByEdgeComponentNumber(mask, edgeToComponent);
        objectives[0] = getFitnessOfComponent(mask, edgeToComponent, x);

        // objectives[0] = info[0];
        // int[] sizes = getComponentSize(edgeToComponent, (int) info[1]);
        // System.out.println("vertices = " + sizes[0]);
        // objectives[1] = sizes[0];
        // objectives[1] = sizes[1];

        return objectives;
    }
    
    public List<Boolean> getComponentByNumber(List<Boolean> mask, int component) {
        List<Integer> edgeToComponent = findComponents(mask);

        List<Boolean> componentMask = new ArrayList<>(mask.size());
        for (int i = 0; i < edges.size(); ++i) {
            if (edgeToComponent.get(i) == component) {
                componentMask.add(true);
            } else {
                componentMask.add(false);
            }
        }

        return componentMask;
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

    // returns: [heaviestComponentFitness, heaviestComponentNumber]
    public double[] getHeaviestComponentInfo(List<Boolean> mask) {
        List<Integer> edgeToComponent = findComponents(mask);
        return getHeaviestComponentInfo(mask, edgeToComponent);
    }

    // returns: [heaviestComponentFitness, heaviestComponentNumber]
    private double[] getHeaviestComponentInfo(List<Boolean> mask,
                                              List<Integer> edgeToComponent) {
        double[] componentInfo = new double[2];

        double fitness = Double.NEGATIVE_INFINITY;
        double component = 0.0;
        for (int i = 1; i <= componentsNumber; ++i) {
            double newFitness = getFitnessOfComponent(mask, edgeToComponent, i);
            if (fitness < newFitness) {
                fitness = newFitness;
                component = i;
            }
        }

        componentInfo[0] = fitness;
        componentInfo[1] = component;

        return componentInfo;
    }

    // returns: [numberOfEdges, numberOfVertices]
    private int[] getComponentSize(List<Integer> edgeToComponent, int component) {
        int[] componentSize = new int[2];
        Map<V, Boolean> isCountedVertex = new HashMap<>();
        for (int i = 0; i < edges.size(); ++i) {
            if (edgeToComponent.get(i) == component) {
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

    private List<Integer> findComponents() {
        List<Boolean> wholeGraphMask = new ArrayList<>(Collections.nCopies(edges.size(), true));
        return findComponents(wholeGraphMask);
    }

    // Use this field carefully as it updates each time with calling 'findComponents(...)'.
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

    // May be called only from 'findComponents(...)' function.
    private void dfs(int edgeNumber, int componentNumber, List<Boolean> mask, List<Integer> edgeToComponent) {
        // todo get rid of double checking that is optimize operations
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

        V second = edge.second;
        for (Map.Entry<V, S> entry : edgesAsMap.get(second).entrySet()) {
            V secondTo = entry.getKey();
            int indexOfEdge = edgesToIndex.get(second).get(secondTo);
            // if edge in graph and isn't counted then
            if (mask.get(indexOfEdge) && edgeToComponent.get(indexOfEdge) != componentNumber) {
                dfs(indexOfEdge, componentNumber, mask, edgeToComponent);
            }
        }
    }

    private double getFitnessOfComponent(List<Boolean> mask, List<Integer> edgeToComponent, int component) {
        double fitness = 0.0;
        Set<S> uniqueSignals = new HashSet<>();
        for (int i = 0; i < edges.size(); ++i) {
            if (mask.get(i) && edgeToComponent.get(i) == component) {
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

    public boolean isNewEdgeInConnectedComponent(List<Boolean> mask, int edgeNumber) {
        Edge edge = edges.get(edgeNumber);
        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                Edge edgeCompareTo = edges.get(i);
                if (edge.first.equals(edgeCompareTo.first) || edge.first.equals(edgeCompareTo.second) ||
                    edge.second.equals(edgeCompareTo.first) || edge.second.equals(edgeCompareTo.second)) {
                    return true;
                }
            }
        }

        return false;
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

    // Returns array of the form:
    // <order_of_connected_component_by_size> -> <number_of_connected_component>
    private int[] orderToNumber() {
        // componentNumberToSize:
        List<ComponentInfo> components = new ArrayList<>(componentsNumberInWholeGraph);
        for (int i = 0; i < componentsNumberInWholeGraph; ++i) {
            components.add(i, new ComponentInfo(i, 0));
        }

        for (int i = 0; i < edges.size(); ++i) {
            int componentNumber = edgeToComponentInWholeGraph.get(i);
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

    /**
     * FOR PRINTING OUTPUT
     */

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
            result.add(new Pair<>(vertex, verticesToSignals.get(vertex)));
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

}
