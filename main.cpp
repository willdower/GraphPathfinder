//
// Created by William on 10/05/2020.
//

#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <limits>
#include <algorithm>
#include <array>
#include <unordered_map>
#include <unordered_set>

struct Edge {
    int source;
    int dest;
    double weight;

    Edge(int pSource, int pDest, double pWeight) {
        source = pSource;
        dest = pDest;
        weight = pWeight;
    }
};

class Graph {
public:
    std::unordered_set<Edge> incidentEdges(int vertex) {
        return adjacencySet[vertex];
    }

    void insertEdge(int source, int dest, double weight) {
        adjacencySet[source].insert(Edge(source, dest, weight));
    }

private:
    std::unordered_map<int, std::unordered_set<Edge>> adjacencySet;
};

int main (int argc, char **argv) {



    return 0;
}