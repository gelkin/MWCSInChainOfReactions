package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

import ec.BreedingPipeline;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Parameter;
import ec.vector.BitVectorIndividual;
import ec.vector.BitVectorSpecies;
import ec.vector.VectorDefaults;

import java.util.ArrayList;
import java.util.List;

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

            List<Boolean> indAsList = new ArrayList(i.genome.length);
            for (int k = 0; k < i.genome.length; ++k) {
                indAsList.add(i.genome[k]);
            }

            for(int x = 0; x < i.genome.length; ++x) {
                /*if (state.random[thread].nextBoolean(species.mutationProbability(x))) {

                    if (!i.genome[x]) {
                        if (graph.isNewEdgeInConnectedComponent(indAsList, x)) {
                            i.genome[x] = true;
                        }
                    } else {
                        i.genome[x] = false;
                    }
                }*/
                // Ordinary mutation
                if (state.random[thread].nextBoolean(species.mutationProbability(x))) {
                    i.genome[x] = !i.genome[x];
                }
            }
            // it's a "new" individual, so it's no longer been evaluated
            i.evaluated = false;
        }

        return n;
    }

}