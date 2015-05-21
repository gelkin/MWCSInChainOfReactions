package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import com.sun.org.apache.xpath.internal.operations.Bool;

import java.io.PrintWriter;
import java.util.*;

/**
 *
 * @param <V> - vertex type
 * @param <S> - signal type
 */
public class Graph<V extends Comparable<V>, S> implements Cloneable {
    private final LinkedHashMap<V, LinkedHashMap<V, S> > edgesAsMap;
    public final LinkedHashMap<V, S> verticesToSignals;
    public final Map<S, Double> signals;
    public final List<Edge<V>> edges;
    public final Map<V, Map<V, Integer>> edgesToIndex;
    private final List<Integer> edgeToComponentInWholeGraph;
    private final int componentsNumberInWholeGraph;
    private int[] orderToNumberForComponents;
    private Map<V, Integer> verticesToDegrees;
    private double maxDegree;

    private static final double SPECIAL_VERTEX_WEIGHT = -1.38033632932348;

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
        initVerticesToDegrees();
        maxDegree = (double) Collections.max(verticesToDegrees.values());
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

    private static final int NUM_OBJECTIVES = 2;
    private double[] multiObjectiveFitness(List<Boolean> mask,
                                           List<Integer> edgeToComponent) {
        double[] objectives = new double[NUM_OBJECTIVES];

        // #1
        // int x = getBiggestByEdgeComponentNumber(mask, edgeToComponent);
        // objectives[0] = getFitnessOfComponent(mask, edgeToComponent, x);

        // #2
        // double[] info = getHeaviestComponentInfo(mask, edgeToComponent);
        // objectives[0] = info[0];
        // int[] sizes = getComponentSize(edgeToComponent, (int) info[1]);
        // System.out.println("vertices = " + sizes[0]);
        // objectives[1] = sizes[0];
        // objectives[2] = sizes[1];

        // #3 in whole graph
        objectives[0] = componentsNumber;
        objectives[1] = getFitnessOfComponent(mask, edgeToComponent, 0);

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

    public int getBiggestByEdgeComponentNumber(List<Boolean> mask) {
        return getBiggestByEdgeComponentNumber(mask, findComponents(mask));
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

    public List<Integer> safeFindComponents(boolean[] mask){
        List<Boolean> maskAsList = toList(mask);

        return safeFindComponents(maskAsList);
    }

    public List<Integer> safeFindComponents(List<Boolean> mask) {
        int tmp = componentsNumber;
        List<Integer> edgeToComponent = findComponents(mask);
        componentsNumber = tmp;

        return edgeToComponent;
    }

    private List<Boolean> toList(boolean[] mask) {
        List<Boolean> maskAsList = new ArrayList(mask.length);
        for (int k = 0; k < mask.length; ++k) {
            maskAsList.add(mask[k]);
        }
        return maskAsList;
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

    // if (component == 0) return fitness of whole subgraph
    public double getFitnessOfComponent(List<Boolean> mask, List<Integer> edgeToComponent, int component) {

        double fitness = 0.0;
        Set<S> uniqueSignals = new HashSet<>();
        for (int i = 0; i < edges.size(); ++i) {
            if (mask.get(i) && ((component == 0) || edgeToComponent.get(i) == component)) {
                Edge<V> edge = edges.get(i);
                uniqueSignals.add(edgesAsMap.get(edge.first).get(edge.second));

                // TODO additional condition: count once only positive weights
                /*if (signals.get(verticesToSignals.get(edge.first)) < 0) {
                    fitness += signals.get(verticesToSignals.get(edge.first));
                } else*/ if (signals.get(verticesToSignals.get(edge.first)) != SPECIAL_VERTEX_WEIGHT) {
                    uniqueSignals.add(verticesToSignals.get(edge.first));
                } else {
                    fitness += SPECIAL_VERTEX_WEIGHT;
                }

                // TODO additional condition: count once only positive weights
                /*if (signals.get(verticesToSignals.get(edge.second)) < 0) {
                    fitness += signals.get(verticesToSignals.get(edge.second));
                } else*/ if (signals.get(verticesToSignals.get(edge.second)) != SPECIAL_VERTEX_WEIGHT) {
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

    public boolean isNewEdgeInConnectedComponent(boolean[] mask, int edgeNumber) {
        List<Boolean> maskAsList = toList(mask);

        Edge<V> edge = edges.get(edgeNumber);

        V first = edge.first;
        for (Map.Entry<V, S> entry : edgesAsMap.get(first).entrySet()) {
            V firstTo = entry.getKey();
            int indexOfEdge = edgesToIndex.get(first).get(firstTo);
            if (maskAsList.get(indexOfEdge)){
                return false;
            }
        }

        V second = edge.second;
        for (Map.Entry<V, S> entry : edgesAsMap.get(second).entrySet()) {
            V secondTo = entry.getKey();
            int indexOfEdge = edgesToIndex.get(second).get(secondTo);
            if (maskAsList.get(indexOfEdge)) {
                return false;
            }
        }

        return true;
    }

    // Returns sum of 'x' vertices degrees, which are not in mask
    public double getEdgeDegree(boolean[] mask, int x) {
        Edge edge = edges.get(x);
        int resDegree = 0;

        boolean isVertexNeeded = true;
        for (Map.Entry<V, S> entry : edgesAsMap.get(edge.first).entrySet()) {
            if (mask[edgesToIndex.get(edge.first).get(entry.getKey())]) {
                isVertexNeeded = false;
                break; // edge.first is already in graph
            }
        }

        if (isVertexNeeded) {
            resDegree += verticesToDegrees.get(edge.first);
        }

        isVertexNeeded = true;
        for (Map.Entry<V, S> entry : edgesAsMap.get(edge.second).entrySet()) {
            if (mask[edgesToIndex.get(edge.second).get(entry.getKey())]) {
                isVertexNeeded = false;
                break; // edge.first is already in graph
            }
        }

        if (isVertexNeeded) {
            resDegree += verticesToDegrees.get(edge.second);
        }

        return ((double) resDegree);
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

    // Returns map of <vertex> -> <degree>
    private void initVerticesToDegrees() {
        verticesToDegrees = new HashMap<>(verticesToSignals.size());
        for (V vertex : verticesToSignals.keySet()) {
            verticesToDegrees.put(vertex, edgesAsMap.get(vertex).size());
        }
    }

    public double getMaxDegree() {
        return maxDegree;
    }

    public List<V> getComponentVertices(int componentNumber, List<Integer> edgeToComponent) {
        Set<V> vertices = new HashSet<>();
        for (int i = 0; i < edges.size(); ++i) {
            if (edgeToComponent.get(i) == componentNumber) {
                Edge<V> edge = edges.get(i);
                vertices.add(edge.first);
                vertices.add(edge.second);
            }
        }

        List<V> verticesAsList = new ArrayList<>(vertices.size());
        verticesAsList.addAll(vertices);
        return verticesAsList;
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

    // Get info about subgraph presented by mask:
    public Map<V, S> getVerticesToSignalsWithNaN(List<Boolean> mask) {
        Map<V, S> verticesToSignalsWithNaN = new LinkedHashMap<>(verticesToSignals.size());
        for (Map.Entry<V, S> entry : verticesToSignals.entrySet()) {
            verticesToSignalsWithNaN.put(entry.getKey(), (S) "NaN");
        }

        Set<V> vertices = new HashSet<>();
        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                vertices.add(edges.get(i).first);
                vertices.add(edges.get(i).second);
            }
        }

        Iterator<V> it = vertices.iterator();
        while (it.hasNext()) {
            V vertex = it.next();
            verticesToSignalsWithNaN.put(vertex, verticesToSignals.get(vertex));
        }

        return verticesToSignalsWithNaN;
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

    public List<Pair<Edge<V>, S>> getEdgesToSignalsWithNaN(List<Boolean> mask) {
        List<Pair<Edge<V>, S>> result = new ArrayList<>();

        for (int i = 0; i < mask.size(); ++i) {
            Edge e = edges.get(i);
            if (mask.get(i)) {
                result.add(new Pair(e, edgesAsMap.get(e.first).get(e.second)));
            } else {
                result.add(new Pair(e, "NaN"));
            }
        }

        return result;
    }

    /**
      Dijkstra
      */

    // Argument 'parent' is assigned and 'returned'
    public Map<V, Double> dijkstra(V start, Map<S, Boolean> contextMap, Map<V, V> parent) {
        Map<V, Double> shortestPath = new HashMap<>(verticesToSignals.keySet().size());
        for (V vertex : verticesToSignals.keySet()) {
            shortestPath.put(vertex, -Double.MAX_VALUE);
        }
        shortestPath.put(start, 0.0);

        PriorityQueue<Map.Entry<V, Double>> queue = new PriorityQueue<>(new Comparator<Map.Entry<V, Double>>() {
            @Override
            public int compare(Map.Entry o1, Map.Entry o2) {
                Double d1 = (Double) o1.getValue();
                Double d2 = (Double) o2.getValue();
                // big double -> come first
                if (d1 < d2) {
                    return 1;
                } else if (d1 > d2) {
                    return -1;
                } else {
                    return ((Comparable<V>) o1.getKey()).compareTo((V) o2.getKey());
                }
            }
        });
        for (Map.Entry<V, Double> entry : shortestPath.entrySet()) {
            queue.add(entry);
        }

        Map<V, Boolean> reached = new HashMap<>(verticesToSignals.size());
        for (V v : verticesToSignals.keySet()) {
            reached.put(v, false);
        }

        //todo
        int counter = 0;
        Map.Entry<V, Double> minVertex;
        while ((minVertex = queue.poll()) != null) {
            System.out.println(start + " !!! " + shortestPath.get(start));
            System.out.println(minVertex.getKey() + " <-> " + minVertex.getValue());
            if (minVertex.getValue() == -Double.MAX_VALUE) {
                // todo
                System.out.println("Met not connected part");
                break;
            }

            // todo
            System.out.println(minVertex.getKey());

            if (!reached.get(minVertex.getKey())) {
                reached.put(minVertex.getKey(), true);

                // TODO! pay attention to all relatives of vertex in order not to count some edge or edge twice
                for (final Map.Entry<V, S> toEdge : edgesAsMap.get(minVertex.getKey()).entrySet()) {
                    double additionalPathWeight = 0.0;
                    S edgesSignal = edgesAsMap.get(minVertex.getKey()).get(toEdge.getKey());
                    if (signals.get(edgesSignal) < 0.0 && !contextMap.get(edgesSignal)) {
                        additionalPathWeight += signals.get(edgesSignal);
                    }

                    if (signals.get(verticesToSignals.get(toEdge.getKey())) < 0.0
                            && !contextMap.get(verticesToSignals.get(toEdge.getKey()))) {
                        additionalPathWeight += signals.get(verticesToSignals.get(toEdge.getKey()));
                    }

                    double newPathWeight = additionalPathWeight + shortestPath.get(minVertex.getKey());
                    if (newPathWeight > shortestPath.get(toEdge.getKey())) {
                        shortestPath.put(toEdge.getKey(), newPathWeight);
                        // Add it again, as add_&&_poll is faster then 'remove()'
                        queue.add(new Map.Entry<V, Double>() {
                            @Override
                            public V getKey() {
                                return toEdge.getKey();
                            }

                            @Override
                            public Double getValue() {
                                return signals.get(toEdge.getValue());
                            }

                            @Override
                            public Double setValue(Double value) {
                                return value;
                            }
                        });
                    }

                    parent.put(minVertex.getKey(), toEdge.getKey());
                }
            }
        }

        return shortestPath;
    }

    public Map<S, Boolean> getCountedSignalAsMap(List<Boolean> mask) {
        Map<S, Boolean> maskAsMap = new HashMap<>(signals.size());
        for (S key : signals.keySet()) {
            maskAsMap.put(key, false);
        }

        for (int i = 0; i < mask.size(); ++i) {
            if (mask.get(i)) {
                Edge edge = edges.get(i);
                maskAsMap.put(edgesAsMap.get(edge.first).get(edge.second), true);
                maskAsMap.put(verticesToSignals.get(edge.first), true);
                maskAsMap.put(verticesToSignals.get(edge.second), true);
            }
        }

        return maskAsMap;
    }

}
