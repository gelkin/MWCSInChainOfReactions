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

        BitVectorIndividual ind = (BitVectorIndividual) individual;

        List<Boolean> indAsList = graph.toList(ind.genome);

        // Multi
        if (!(ind.fitness instanceof MultiObjectiveFitness)) {
            state.output.fatal("Whoa!  It's not a MultiObjectiveFitness!!!", null);
        }

        double[] objectives = graph.fitness(indAsList, 0);
        ((MultiObjectiveFitness) ind.fitness).setObjectives(state, objectives);

        ind.evaluated = true;
    }

}
