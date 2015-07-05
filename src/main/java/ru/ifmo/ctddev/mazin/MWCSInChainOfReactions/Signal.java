package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions;

public class Signal {
    private double value;

    public Signal(double value) {
        this.value = value;
    }

    public double getValue() {
        return value;
    }

    @Override
    public int hashCode() {
        return (new Double(value)).hashCode();
    }

    @Override
    public boolean equals(Object other) {
        return (hashCode() == other.hashCode());
    }

    @Override
    public String toString() {
        return Integer.toHexString(super.toString().hashCode());
    }
}
