package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import com.sun.javafx.collections.MappingChange;
import ec.BreedingPipeline;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Parameter;
import ec.vector.BitVectorIndividual;
import ec.vector.BitVectorSpecies;
import ec.vector.VectorDefaults;

import java.util.*;

public class GraphMutatorPipeline extends BreedingPipeline {
    public static final String P_MYMUTATION = "my-mutation";
    private static Graph graph;

    private static boolean isGraphSet = false;
    public static void setGraph(Graph graph) {
        if (!isGraphSet) {
            GraphMutatorPipeline.graph = graph;
            isGraphSet = true;
        }
    }

    // We have to specify a default base, even though we never use it
    public Parameter defaultBase() {
        return VectorDefaults.base().push(P_MYMUTATION);
    }

    public static final int NUM_SOURCES = 1;

    // Return 1 -- we only use one source
    public int numSources() {
        return NUM_SOURCES;
    }

    // We're supposed to create a most _max_ and at least _min_ individuals,
    // drawn from our source and mutated, and stick them into slots in inds[]
    // starting with the slot inds[start].  Let's do this by telling our
    // source to stick those individuals into inds[] and then mutating them
    // right there.
    public int produce(final int min,
                       final int max,
                       final int start,
                       final int subpopulation,
                       final Individual[] inds,
                       final EvolutionState state,
                       final int thread)
    {
        // grab individuals from our source and stick 'em right into inds.
        // we'll modify them from there
        int n = sources[0].produce(min,max,start,subpopulation,inds,state,thread);


        // should we bother?
        if (!state.random[thread].nextBoolean(likelihood)) {
            return reproduce(n, start, subpopulation, inds, state, thread, false);  // DON'T produce children from source -- we already did
        }

        // clone the individuals if necessary -- if our source is a BreedingPipeline
        // they've already been cloned, but if the source is a SelectionMethod, the
        // individuals are actual individuals from the previous population
        if (!(sources[0] instanceof BreedingPipeline)) {
            for (int q = start; q < n + start; ++q) {
                inds[q] = (Individual) (inds[q].clone());
            }
        }

        // Check to make sure that the individuals are BitVectorIndividuals and
        // grab their species.  For efficiency's sake, we assume that all the
        // individuals in inds[] are the same type of individual and that they all
        // share the same common species -- this is a safe assumption because they're
        // all breeding from the same subpopulation.

        if (!(inds[start] instanceof BitVectorIndividual)) { // uh oh, wrong kind of individual
            state.output.fatal("OurMutatorPipeline didn't get an BitVectorIndividual." +
                    "The offending individual is in subpopulation " + subpopulation + " and it's:" + inds[start]);
        }
        BitVectorSpecies species = (BitVectorSpecies) (inds[start].species);

        // mutate 'em!
        for(int q = start; q < (n + start); ++q) {
            BitVectorIndividual i = (BitVectorIndividual)inds[q];

            // first part
            // TODO connect all connected components with some possibility
            // Let's connect two random connected component
            List<Integer> edgesToComponents = graph.safeFindComponents(i.genome);
            List<Integer> componentNumbers = getLeftComponentsNumbers(edgesToComponents);
            if (componentNumbers.size() > 1) {
                int firstComp = state.random[thread].nextInt(componentNumbers.size());
                int secondComp;
                if ((secondComp = state.random[thread].nextInt(componentNumbers.size() - 1)) >= firstComp) {
                    ++secondComp; // so that secondComp != firstComp
                }
                List<String> firstVertices = graph.getComponentVertices(firstComp, edgesToComponents);
                List<String> secondVertices = graph.getComponentVertices(secondComp, edgesToComponents);
                // It doesn't matter which vertex to chose as start-vertex, because all vertices and edges
                // inside one component are counted as zero-weight.
                // So lets take zero vertex:
                List<Integer> pathEdges = getPathEdges(firstVertices.get(0), secondVertices, i.genome);

                // lets add new edges to individual
                for (int j = 0; j < pathEdges.size(); ++j) {
                    i.genome[pathEdges.get(j)] = true;
                }
            }

            // second part
            for(int x = 0; x < i.genome.length; ++x) {
                // #1 Add new edge which is not in any connected component with lower probability
                if (!i.genome[x]) {
                    double prob = species.mutationProbability(x);
                    if (graph.isNewEdgeInConnectedComponent(i.genome, x)) {
                        prob *= prob;
                    }

                    if (state.random[thread].nextBoolean(prob)) {
                        i.genome[x] = true;
                    }
                } else {
                    if (state.random[thread].nextBoolean(species.mutationProbability(x))) {
                        i.genome[x] = false;
                    }
                }

                // #2 Connection && degree
                /*
                add private static final double DEGREE_COEF = 5.0; // TODO must be (< 1.0)
                double edgeDegree = graph.getEdgeDegree(i.genome, x) / graph.getMaxDegree();
                if (!i.genome[x]) {
                    double prob = species.mutationProbability(x);
                    if (graph.isNewEdgeInConnectedComponent(i.genome, x)) {
                        prob *= prob;
                    }

                    if (state.random[thread].nextBoolean(prob * (1.0 + DEGREE_COEF * edgeDegree))) {
                        i.genome[x] = true;
                    }
                } else {
                    if (state.random[thread].nextBoolean(species.mutationProbability(x))) {
                        i.genome[x] = false;
                    }
                }*/

                // #3 Ordinary mutation
                /*
                if (state.random[thread].nextBoolean(species.mutationProbability(x))) {
                    i.genome[x] = !i.genome[x];
                }
                */

            }

            // it's a "new" individual, so it's no longer been evaluated
            i.evaluated = false;
        }

        return n;
    }

    private List<Integer> getLeftComponentsNumbers(List<Integer> ind) {
        Set<Integer> numbers = new HashSet<>();
        for (int i = 0; i < ind.size(); ++i) {
            if (ind.get(i) != 0) {
                numbers.add(ind.get(i));
            }
        }

        List<Integer> numbersAsList = new ArrayList<>(numbers.size());
        numbersAsList.addAll(numbers);

        return numbersAsList;
    }

    private List<Integer> getPathEdges(String start, List<String> secondComponent, boolean[] mask) {
        List<Boolean> maskAsList = new ArrayList<>(mask.length);
        for (int i = 0; i < mask.length; ++i) {
            maskAsList.add(mask[i]);
        }

        Map<String, String> parent = new HashMap<>();
        Map<String, Double> shortestPath = graph.dijkstra(start, graph.getCountedSignalAsMap(maskAsList), parent);

        String maxVertex = start;
        double maxValue = -Double.MAX_VALUE;
        for (String v : secondComponent) {
            if (shortestPath.get(v) > maxValue) {
                maxValue = shortestPath.get(v);
                maxVertex = v;
            }
        }
        // components cannot be connected
        if (maxVertex.equals(start)) {
            return new ArrayList<>();
        }

        List<Edge> pathEdges = new ArrayList<>();
        String parentVertex = maxVertex;
        while (!parentVertex.equals(start)) {
            pathEdges.add(new Edge(parent.get(parentVertex), parentVertex));
            parentVertex = parent.get(parentVertex);
        }

        List<Integer> pathEdgesNumbers = new ArrayList<>();
        for (Edge e : pathEdges) {
            pathEdgesNumbers.add((Integer) ((Map) graph.edgesToIndex.get(e.first)).get(e.second));
        }

        return pathEdgesNumbers;
    }

}