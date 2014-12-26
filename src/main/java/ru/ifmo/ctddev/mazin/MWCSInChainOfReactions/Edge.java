package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

public class Edge<V> {
    V first;
    V second;

    public Edge(V first, V second) {
        this.first = first;
        this.second = second;
    }
}
