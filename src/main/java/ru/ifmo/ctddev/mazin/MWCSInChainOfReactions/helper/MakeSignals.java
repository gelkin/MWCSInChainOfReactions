package ru.ifmo.ctddev.mazin.MWCSInChainOfReactions.helper;

import ru.ifmo.ctddev.mazin.MWCSInChainOfReactions.*;

import java.io.*;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;


public class MakeSignals {
    private static final String VERTICES_FILE = "./src/main/resources/old_resources/old_nodes.txt";
    private static final String EDGES_FILE = "./src/main/resources/old_resources/old_edges.txt";

    public static final String NODES_PATH = "./src/main/resources/nodes.txt";
    public static final String EDGES_PATH = "./src/main/resources/edges.txt";
    public static final String SIGNALS_PATH = "./src/main/resources/signals.txt";

    public static void main(String[] args) {

        PrintWriter nodesWriter = null;
        PrintWriter edgesWriter = null;
        PrintWriter signalsWriter = null;
        try {
            nodesWriter = new PrintWriter(NODES_PATH, "UTF-8");
            edgesWriter = new PrintWriter(EDGES_PATH, "UTF-8");
            signalsWriter = new PrintWriter(SIGNALS_PATH, "UTF-8");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        }

        // Vertices first:
        Map<Signal, Double> signals = new HashMap<>();
        LinkedHashMap<String, Signal> verticesToSignal = new LinkedHashMap<>();
        try(BufferedReader br = new BufferedReader(new FileReader(VERTICES_FILE))) {
            String line = br.readLine();

            while (line != null) {
                String[] vertexInfo = line.split("\\t");
                String vertex = vertexInfo[0];
                double weight = Double.parseDouble(vertexInfo[1]);

                Signal signal = new Signal(weight);
                // should replace same signal's value
                signals.put(signal, weight);

                verticesToSignal.put(vertex, signal);

                line = br.readLine();
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Edges
        LinkedHashMap<String, LinkedHashMap<String, Signal>> edgesAsMap = new LinkedHashMap<>();
        List<Edge<String>> edges = new ArrayList<>();
        Map<String, Map<String, Integer>> edgesToIndex = new HashMap<>();

        // todo check if I add one signal key only once
        int edgeNumber = 0;
        try(BufferedReader br = new BufferedReader(new FileReader(EDGES_FILE))) {
            String line = br.readLine();

            while (line != null) {
                String[] edgeInfo = line.split("\\t");
                String from = edgeInfo[0];
                String to = edgeInfo[1];

                double weight = Double.parseDouble(edgeInfo[2]);
                Signal signal = new Signal(weight);
                signals.put(signal, weight);

                if (!edgesAsMap.containsKey(from)) {
                    edgesAsMap.put(from, new LinkedHashMap<String, Signal>());
                    edgesToIndex.put(from, new HashMap<String, Integer>());
                }

                if (!edgesAsMap.containsKey(to)) {
                    edgesAsMap.put(to, new LinkedHashMap<String, Signal>());
                    edgesToIndex.put(to, new HashMap<String, Integer>());
                }

                // if this edge wasn't added
                if (!edgesAsMap.get(from).containsKey(to) && !edgesAsMap.get(to).containsKey(from)) {
                    edgesAsMap.get(from).put(to, signal);
                    edgesAsMap.get(to).put(from, signal);

                    edges.add(new Edge<>(from, to));

                    edgesToIndex.get(from).put(to, edgeNumber);
                    edgesToIndex.get(to).put(from, edgeNumber);

                    ++edgeNumber;
                }

                line = br.readLine();
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Make new resources:

        Map<Signal, Boolean> isAddedSignal = new HashMap<>();
        // write vertices
        for (Map.Entry<String, Signal> vertex : verticesToSignal.entrySet()) {
            String signalName = vertex.getValue().toString();
            nodesWriter.write(vertex.getKey()+ "\t" + signalName + System.lineSeparator());

            if ("705d88c0".equals(signalName)) {
                int x = 0xdeadbeef;
            }

            if (!isAddedSignal.containsKey(vertex.getValue())) {

                System.out.println(signalName + "\t" + vertex.getValue());

                signalsWriter.write(signalName + "\t" + signals.get(vertex.getValue()) + System.lineSeparator());
                isAddedSignal.put(vertex.getValue(), true);
            }
        }

        // write edges
        Map<String, Set<String>> areEdgesWritten = new HashMap<>();
        for (Map.Entry<String, LinkedHashMap<String, Signal>> edge : edgesAsMap.entrySet()) {
            areEdgesWritten.put(edge.getKey(), new HashSet<String>());
            for (Map.Entry<String, Signal> edgeTo : edge.getValue().entrySet()) {
                // if we iterated over vertex "edgeTo" then we added {edge, edgeTo}
                if (!areEdgesWritten.containsKey(edgeTo.getKey())) {
                    areEdgesWritten.get(edge.getKey()).add(edgeTo.getKey());
                    String signalName = edgeTo.getValue().toString();
                    edgesWriter.write(edge.getKey() + "\t" +
                                      edgeTo.getKey() + "\t" +
                                      signalName +
                                      System.lineSeparator());

                    // add signal if needed
                    if (!isAddedSignal.containsKey(edgeTo.getValue())) {
                        signalsWriter.write(signalName + "\t" + signals.get(edgeTo.getValue()) + System.lineSeparator());
                        isAddedSignal.put(edgeTo.getValue(), true);
                    }
                }
            }
        }

        nodesWriter.flush();
        edgesWriter.flush();
        signalsWriter.flush();
    }

    public static String padRight(String s, int n) {
        return String.format("%1$-" + n + "s", s);
    }

    public static String padLeft(String s, int n) {
        return String.format("%1$" + n + "s", s);
    }

    public static String getHash(String str) {
        try {
            MessageDigest messageDigest = MessageDigest.getInstance("MD5");
            messageDigest.update(str.getBytes());
            return new String(messageDigest.digest());
        } catch (NoSuchAlgorithmException e) {
            e.printStackTrace();
            return null;
        }

    }
}
