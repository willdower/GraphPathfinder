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
    std::unordered_set<Edge*> incidentEdges(int vertex) {
        return adjacencySet[vertex];
    }

    void insertEdge(int source, int dest, double weight) {
        adjacencySet[source].insert(new Edge(source, dest, weight));
    }

    void insertEdge(Edge *newEdge) {
        adjacencySet[newEdge->source].insert(newEdge);
    }

    void removeEdge(Edge *edgeToRemove) {
        adjacencySet[edgeToRemove->source].erase(edgeToRemove);
    }

    int getNumberOfVertices() {
        return vertices;
    }

    int getNumberOfEdges() {
        return edges;
    }

    std::vector<int> vertexSet() {
        std::vector<int> returnVector;
        for (auto const & vertex : adjacencySet) {
            returnVector.push_back(vertex.first);
        }
        return returnVector;
    }

    Graph(int pVertices, int pEdges) {
        vertices = pVertices;
        edges = pEdges;
    }

private:
    std::unordered_map<int, std::unordered_set<Edge*>> adjacencySet;
    int vertices;
    int edges;
};

class VertexComparator {
public:
    static double *distArr;

    std::priority_queue<int> test;

    bool operator () (int a, int b) {
        test.pop();
        return distArr[a] > distArr[b];
    }
};

class CustomPriorityQueue {
public:
    void heapify() {
        std::make_heap(queueVector.begin(), queueVector.end(), VertexComparator());
    }

    void push(int toPush) {
        std::push_heap(queueVector.begin(), queueVector.end(), VertexComparator());
    }

    void pop() {
        std::pop_heap(queueVector.begin(), queueVector.end());
        queueVector.pop_back();
    }

    bool empty() {
        return queueVector.empty();
    }

    int top() {
        return queueVector.front();
    }

    CustomPriorityQueue() = default;

    explicit CustomPriorityQueue(std::vector<int> const & initialVector) {
        queueVector = initialVector;
        heapify();
    }
private:
    std::vector<int> queueVector;
};

std::vector<int> dijkstra(Graph *graph, int source, int sink) {
    std::vector<int> returnPath;
    auto dist = new double [graph->getNumberOfEdges()];
    for (int i=0;i<graph->getNumberOfEdges();i++) {
        dist[i] = std::numeric_limits<double>::max();
    }
    auto prev = new int [graph->getNumberOfEdges()]{0};

    dist[source] = 0;

    VertexComparator::distArr = dist;
    CustomPriorityQueue nodeQueue(graph->vertexSet());

    while (!nodeQueue.empty()) {
        int current = nodeQueue.top();
        nodeQueue.pop();

        std::unordered_set<Edge*> incidentEdges = graph->incidentEdges(current);
        for (auto edge : incidentEdges) {
            double originalDistance = dist[edge->dest];
            double newDistance = (dist[edge->source] + edge->weight);
            if (originalDistance > newDistance) {
                dist[edge->dest] = dist[edge->source] + edge->weight;
                prev[edge->dest] = edge->source;
                nodeQueue.heapify();
            }
        }
    }

    int current = sink;
    while (current != source) {
        returnPath.push_back(current);
        if (current != sink) {
        }
        current = prev[current];
    }
    returnPath.push_back(current);
    std::reverse(returnPath.begin(), returnPath.end());

    std::cout << dist[sink] << std::endl;

    delete [] dist;
    delete [] prev;

    return returnPath;
}

double printPath(std::vector<int> const & path) {
    for (auto node : path) {
        std::cout << "->" << node;
    }
    std::cout << std::endl;
}

double *VertexComparator::distArr = nullptr;

int main (int argc, char **argv) {

    if (argc != 2) {
        std::cout << "Incorrect number of arguments." << std::endl;
        exit(1);
    }

    int vertices, edges;

    std::ifstream input(argv[1]);

    input >> vertices;
    input >> edges;

    if (vertices == 0 || edges == 0) {
        std::cout << "No edges or no vertices" << std::endl;
        exit(2);
    }

    auto *originalGraph = new Graph(vertices, edges);

    for (int i=0;i<edges;i++) {
        int start, end;
        double weight;
        input >> start;
        input >> end;
        input >> weight;
        originalGraph->insertEdge(start, end, weight);
    }

    int source, destination, k;

    input >> source;
    input >> destination;
    input >> k;

    printPath(dijkstra(originalGraph, source, destination));

    auto *pathArray = new std::vector<Edge*>[vertices];

    return 0;
}