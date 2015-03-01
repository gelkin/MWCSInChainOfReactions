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
    private static Graph graph;

    private static boolean isGraphSet = false;
    public static void setGraph(Graph graph) {
        if (!isGraphSet) {
            MWCGProblem.graph = graph;
            isGraphSet = true;
        }
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


        // Multi
        if (!(ind2.fitness instanceof MultiObjectiveFitness)) {
            state.output.fatal("Whoa!  It's not a MultiObjectiveFitness!!!", null);
        }

        double[] objectives = graph.fitness(ind2AsList, 1);
        ((MultiObjectiveFitness) ind2.fitness).setObjectives(state, objectives);


        /*
        if (!(ind2.fitness instanceof SimpleFitness)) {
            state.output.fatal("Whoa!  It's not a SimpleFitness!!!", null);
        }

        ((SimpleFitness)ind2.fitness).setFitness(state,
                                                 graph.fitness(ind2AsList, 0),
                                                 false);
                                                 */

        ind2.evaluated = true;
    }

}
