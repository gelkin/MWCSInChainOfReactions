package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;


import ec.EvolutionState;
import ec.vector.BitVectorIndividual;
import ec.vector.BitVectorSpecies;
import ec.vector.VectorIndividual;
import ec.vector.VectorSpecies;

import java.util.*;

public class GraphIndividual extends BitVectorIndividual {

    private static Graph graph;

    private static boolean isGraphSet = false;
    public static void setGraph(Graph graph) {
        if (!isGraphSet) {
            GraphIndividual.graph = graph;
            isGraphSet = true;
        }
    }

    @Override
    public void defaultCrossover(EvolutionState state, int thread, VectorIndividual ind) {
        BitVectorSpecies s = (BitVectorSpecies)species;  // where my default info is stored
        BitVectorIndividual i = (BitVectorIndividual) ind;
        boolean tmp;
        int point;

        int len = Math.min(genome.length, i.genome.length);
        if (len != genome.length || len != i.genome.length)
            state.output.warnOnce("Genome lengths are not the same.  Vector crossover will only be done in overlapping region.");



        switch(s.crossoverType)
        {
            // TODO IT IS NOT WHAT IT SEEMS! I don't want to change VertorSpices, so I use
            // TODO 'one'-crossover name for my own purposes.
            case VectorSpecies.C_ONE_POINT:
                boolean[] intersect = new boolean[genome.length];
                for (int j = 0; j < genome.length; ++j) {
                    intersect[j] = genome[j] && i.genome[j];
                }

                // deal with components
                List<Integer> componentsFirst = (List) graph.findComponents(genome).first;
                List<Integer> componentsSecond = (List) graph.findComponents(i.genome).first;

                // put 'genome' and 'i.genome' to be equal to 'intersect'
                System.arraycopy(intersect, 0, genome, 0, genome.length);
                System.arraycopy(intersect, 0, i.genome, 0, i.genome.length);

                // get rid of edges in intersection
                for (int j = 0; j < genome.length; ++j) {
                    if (intersect[j]) {
                        componentsFirst.set(j, 0);
                        componentsSecond.set(j, 0);
                    }
                }

                rewriteIndividualsAccordingToComponents(state, thread, i.genome, componentsFirst);
                rewriteIndividualsAccordingToComponents(state, thread, i.genome, componentsSecond); // 'i.genome' is okay here

                break;
            case VectorSpecies.C_ONE_POINT_NO_NOP:
                point = state.random[thread].nextInt((len / s.chunksize) - 1) + 1;  // so it goes from 1 .. len-1
                for(int x=0;x<point*s.chunksize;x++)
                {
                    tmp = i.genome[x];
                    i.genome[x] = genome[x];
                    genome[x] = tmp;
                }
                break;
            case VectorSpecies.C_TWO_POINT:
            {
//                int point0 = state.random[thread].nextInt((len / s.chunksize)+1);
//                point = state.random[thread].nextInt((len / s.chunksize)+1);
                // we want to go from 0 to len-1
                // so that the only NO-OP crossover possible is point == point0
                // example; len = 4
                // possibilities: a=0 b=0       NOP                             [0123]
                //                                a=0 b=1       swap 0                  [for 1, 2, 3]
                //                                a=0 b=2       swap 0, 1               [for 2, 3]
                //                                a=0 b=3       swap 0, 1, 2    [for 3]
                //                                a=1 b=1       NOP                             [1230]
                //                                a=1 b=2       swap 1                  [for 2, 3, 0]
                //                                a=1 b=3       swap 1, 2               [for 3, 0]
                //                                a=2 b=2       NOP                             [2301]
                //                                a=2 b=3       swap 2                  [for 3, 0, 1]
                //                                a=3 b=3   NOP                         [3012]
                // All intervals: 0, 01, 012, 0123, 1, 12, 123, 1230, 2, 23, 230, 2301, 3, 30, 301, 3012
                point = state.random[thread].nextInt((len / s.chunksize));
                int point0 = state.random[thread].nextInt((len / s.chunksize));
                if (point0 > point) { int p = point0; point0 = point; point = p; }
                for(int x=point0*s.chunksize;x<point*s.chunksize;x++)
                {
                    tmp = i.genome[x];
                    i.genome[x] = genome[x];
                    genome[x] = tmp;
                }
            }
            break;
            case VectorSpecies.C_TWO_POINT_NO_NOP:
            {
                point = state.random[thread].nextInt((len / s.chunksize));
                int point0 = 0;
                do { point0 = state.random[thread].nextInt((len / s.chunksize)); }
                while (point0 == point);  // NOP
                if (point0 > point) { int p = point0; point0 = point; point = p; }
                for(int x=point0*s.chunksize;x<point*s.chunksize;x++)
                {
                    tmp = i.genome[x];
                    i.genome[x] = genome[x];
                    genome[x] = tmp;
                }
            }
            break;
            case BitVectorSpecies.C_ANY_POINT:
                for(int x=0;x<len/s.chunksize;x++)
                    if (state.random[thread].nextBoolean(s.crossoverProbability))
                        for(int y=x*s.chunksize;y<(x+1)*s.chunksize;y++)
                        {
                            tmp = i.genome[y];
                            i.genome[y] = genome[y];
                            genome[y] = tmp;
                        }
                break;
            default:


                state.output.fatal("In valid crossover type in BitVectorIndividual.");
                break;
        }
    }

    private Set<Integer> getLeftComponentsNumbers(List<Integer> ind) {
        Set<Integer> numbers = new HashSet<>();
        for (int i = 0; i < ind.size(); ++i) {
            if (ind.get(i) != 0) {
                numbers.add(ind.get(i));
            }
        }

        return numbers;
    }

    private void rewriteIndividualsAccordingToComponents(EvolutionState state,
                                                         int thread,
                                                         boolean[] indGenome,
                                                         List<Integer> components) {
        Set<Integer> componentsNumbers = getLeftComponentsNumbers(components);

        // (i, b): (b)? put ind3[i] = true: put ind4[i] = true;
        // So every edge from former connected component will be put in only
        // one individual.
        Map<Integer, Boolean> componentToInds = new HashMap<>(componentsNumbers.size());
        for (Integer componentsNumber : componentsNumbers) {
            componentToInds.put(componentsNumber, state.random[thread].nextBoolean());
        }

        // rewrite genome and i.genome
        for (int j = 0; j < genome.length; ++j) {
            if (components.get(j) != 0) {
                int number = components.get(j);
                genome[j] = componentToInds.get(number);
                indGenome[j] = !componentToInds.get(number); // not
            }
        }
    }
}