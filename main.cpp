//
// Created by William on 10/05/2020.
//

#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <chrono>

struct Edge {
    int source;
    int dest;
    double weight;

    bool operator == (const Edge & otherEdge) const {
        return (source == otherEdge.source && dest == otherEdge.dest && weight == otherEdge.weight);
    }

    bool operator != (const Edge & otherEdge) const {
        return !(source == otherEdge.source && dest == otherEdge.dest && weight == otherEdge.weight);
    }

    Edge(int pSource, int pDest, double pWeight) {
        source = pSource;
        dest = pDest;
        weight = pWeight;
    }
};


struct HashFunction {
    size_t operator () (const Edge & k) const {
        return k.dest;
    }
};

class Graph {
public:

    std::unordered_set<Edge, HashFunction> incidentEdges(int vertex) {
        return adjacencySet[vertex];
    }

    void insertEdge(int source, int dest, double weight) {
        adjacencySet[source].insert(Edge(source, dest, weight));
    }

    void removeEdge(Edge edgeToRemove) {
        adjacencySet[edgeToRemove.source].erase(edgeToRemove);
    }

    void removeEdge(int start, int dest) {
        for (auto edge : adjacencySet[start]) {
            if (edge.dest == dest) {
                removeEdge(edge);
                break;
            }
        }
    }

    Edge getEdge(int start, int dest) {
        for (auto edge : adjacencySet[start]) {
            if (edge.dest == dest) {
                return edge;
            }
        }
        return {-1, -1, 0};
    }

    int getNumberOfVertices() {
        return vertices;
    }

    double getPathLength(const std::vector<int> & path, int source) {
        int current = source;
        double totalCost = 0;
        for (auto vertex : path) {
            double weight = getEdge(current, vertex).weight;
            totalCost += weight;
            current = vertex;
        }
        return totalCost;
    }

    void removeLeavingEdges(int vertex) {
        adjacencySet[vertex].clear();
    }

    Graph(int pVertices, int pEdges) {
        vertices = pVertices;
        edges = pEdges;
    }

    explicit Graph(Graph *graph) {
        adjacencySet = graph->adjacencySet;
        vertices = graph->vertices;
        edges = graph->edges;
    }

private:
    std::unordered_map<int, std::unordered_set<Edge, HashFunction>> adjacencySet;
    int vertices;
    int edges;
};

class VertexComparator {
public:
    bool operator () (std::pair<int, double> a, std::pair<int, double> b) {
        return a.second > b.second;
    }
};

class PathComparator {
public:
    bool operator () (std::pair<double, std::vector<int>> const & a, std::pair<double, std::vector<int>> const & b) {
        return a.first < b.first;
    }
};

int *getOptimals(Graph *graph, int source, int sink) {
    source = sink;

    auto dist = new double [graph->getNumberOfVertices()];
    auto visited = new bool [graph->getNumberOfVertices()];
    auto prev = new int [graph->getNumberOfVertices()];

    for (int i=0;i<graph->getNumberOfVertices();i++) {
        dist[i] = std::numeric_limits<double>::max();
        visited[i] = false;
        prev[i] = -1;
    }

    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, VertexComparator> nodeQueue;

    dist[source] = 0;

    std::pair<int, double> currentPair;
    currentPair.first = source;
    currentPair.second = dist[source];
    nodeQueue.push(currentPair);

    while (!nodeQueue.empty()) {
        currentPair = nodeQueue.top();
        nodeQueue.pop();

        visited[currentPair.first] = true;

        for (auto edge : graph->incidentEdges(currentPair.first)) {
            if (!visited[edge.dest] && dist[edge.dest] > dist[edge.source] + edge.weight) {
                dist[edge.dest] = dist[edge.source] + edge.weight;
                prev[edge.dest] = edge.source;
                currentPair.first = edge.dest;
                currentPair.second = dist[edge.dest];
                nodeQueue.push(currentPair);
            }
        }
    }

    delete [] dist;
    delete [] visited;

    return prev;
}

std::vector<int> dijkstra(Graph *graph, int source, int sink, std::unordered_set<int> & visited) {

    auto dist = new double [graph->getNumberOfVertices()];
    auto prev = new int [graph->getNumberOfVertices()];

    for (int i=0;i<graph->getNumberOfVertices();i++) {
        dist[i] = std::numeric_limits<double>::max();
        prev[i] = -1;
    }

    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, VertexComparator> nodeQueue;

    dist[source] = 0;

    std::pair<int, double> currentPair;
    currentPair.first = source;
    currentPair.second = dist[source];
    nodeQueue.push(currentPair);

    auto startTime = std::chrono::system_clock::now();
    while (!nodeQueue.empty()) {
        currentPair = nodeQueue.top();
        nodeQueue.pop();

        visited.insert(currentPair.first);

        for (auto edge : graph->incidentEdges(currentPair.first)) {
            if (!visited.contains(edge.dest) && dist[edge.dest] > dist[edge.source] + edge.weight) {
                dist[edge.dest] = dist[edge.source] + edge.weight;
                prev[edge.dest] = edge.source;
                currentPair.first = edge.dest;
                currentPair.second = dist[edge.dest];
                nodeQueue.push(currentPair);
            }
        }
    }
    auto endTime = std::chrono::system_clock::now();
    std::chrono::duration<double> timeTaken = endTime - startTime;
    std::cout << timeTaken.count() << std::endl;


    int current = sink;
    std::vector<int>returnPath;
    while (current != source) {
        returnPath.push_back(current);
        current = prev[current];
        if (current == -1) {
            break;
        }
    }
    returnPath.push_back(current);
    if (current == -1) {
        returnPath.clear();
    }
    std::reverse(returnPath.begin(), returnPath.end());


    delete [] dist;
    delete [] prev;

    return returnPath;
}

std::vector<int> quickDijkstra(Graph *graph, int source, int sink, std::unordered_set<int> & visited, std::pair<double, std::vector<int>> *optimalPaths, std::unordered_set<int> *optimalPathSet) {

    auto dist = new double [graph->getNumberOfVertices()];
    auto prev = new int [graph->getNumberOfVertices()];

    for (int i=0;i<graph->getNumberOfVertices();i++) {
        dist[i] = std::numeric_limits<double>::max();
        prev[i] = -1;
    }

    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, VertexComparator> nodeQueue;

    dist[source] = 0;

    std::pair<int, double> currentPair;
    currentPair.first = source;
    currentPair.second = dist[source];
    nodeQueue.push(currentPair);

    auto startTime = std::chrono::system_clock::now();
    while (!nodeQueue.empty()) {
        currentPair = nodeQueue.top();
        nodeQueue.pop();

        visited.insert(currentPair.first);
        for (auto edge : graph->incidentEdges(currentPair.first)) {
            if (!visited.contains(edge.dest) && dist[edge.dest] > dist[edge.source] + edge.weight) {
                bool check = true;
                for (auto vertex : visited) {
                    if (optimalPathSet[edge.dest].contains(vertex)) {
                        check = false;
                        break;
                    }
                }
                if (!check) {
                    // Optimal path from this node to the sink contains visited or removed nodes
                    // Regular Dijkstra behaviour
                    dist[edge.dest] = dist[edge.source] + edge.weight;
                    prev[edge.dest] = edge.source;
                    currentPair.first = edge.dest;
                    currentPair.second = dist[edge.dest];
                    nodeQueue.push(currentPair);
                }
                else {
                    // Optimal path is a successful one
                    // Clears all leaving edges and only puts optimal one to destination
                    dist[edge.dest] = dist[edge.source] + edge.weight;
                    prev[edge.dest] = edge.source;
                    graph->removeLeavingEdges(edge.dest);
                    graph->insertEdge(edge.dest, sink, optimalPaths[edge.dest].first);
                    currentPair.first = edge.dest;
                    currentPair.second = optimalPaths[edge.dest].first;
                    if (currentPair.first != sink) {
                        nodeQueue.push(currentPair);
                    }
                }
            }
        }
    }
    auto endTime = std::chrono::system_clock::now();
    std::chrono::duration<double> timeTaken = endTime - startTime;
    std::cout << timeTaken.count() << std::endl;

    int current = sink;
    std::vector<int>returnPath;
    while (current != source) {
        returnPath.push_back(current);
        current = prev[current];
        if (current == -1) {
            break;
        }
    }
    returnPath.push_back(current);
    if (current == -1) {
        returnPath.clear();
    }
    auto it = returnPath.begin() + 1;
    *it = *it * -1; // Flip the node that has been quickened to negative so that information can be passed out of the function
    std::reverse(returnPath.begin(), returnPath.end());


    delete [] dist;
    delete [] prev;

    return returnPath;
}

std::vector<std::pair<double, std::vector<int>>> newYen(Graph *graph, Graph *reverseGraph, int source, int sink, int k) {
    int *optimalPrevs = getOptimals(reverseGraph, source, sink);

    auto *optimalPaths = new std::pair<double, std::vector<int>>[graph->getNumberOfVertices()];
    auto *optimalPathSets = new std::unordered_set<int>[graph->getNumberOfVertices()];
    for (int i=0;i<graph->getNumberOfVertices();i++) {
        int current = i;
        while (current != sink) {
            optimalPaths[i].second.push_back(current);
            optimalPathSets[i].insert(current);
            current = optimalPrevs[current];
            if (current == -1) {
                break;
            }
        }
        optimalPaths[i].second.push_back(current);
        optimalPathSets[i].insert(current);
        optimalPaths[i].first = graph->getPathLength(optimalPaths[i].second, i);
    }

    std::vector<std::pair<double, std::vector<int>>> validPaths;
    std::vector<std::pair<double, std::vector<int>>> candidatePaths;

    std::unordered_set<int> visited;

    std::vector<int> optimalPath = dijkstra(graph, source, sink, visited);

    std::pair<double, std::vector<int>> optimalPathPair;
    optimalPathPair.first = graph->getPathLength(optimalPath, source);
    optimalPathPair.second = optimalPath;
    validPaths.push_back(optimalPathPair);

    auto editableGraph = new Graph(graph);
    for (int i=1;i<=k;i++) {
        for (int j=0;j<=validPaths[i-1].second.size() - 2;j++) {
            visited.clear();
            int spurNode = validPaths[i-1].second[j];
            std::vector<int> rootPath = validPaths[i-1].second;
            auto pos = std::find(rootPath.begin(), rootPath.end(), spurNode);
            rootPath.erase(pos+1, rootPath.end());

            for (auto path : validPaths) {
                auto newPos = std::find(path.second.begin(), path.second.end(), spurNode);
                path.second.erase(newPos+1, path.second.end());
                if (rootPath == path.second) {
                    int start = *newPos;
                    int dest = *(newPos+1);
                    editableGraph->removeEdge(start, dest);
                }
            }

            auto spurPath = quickDijkstra(editableGraph, spurNode, sink, visited, optimalPaths, optimalPathSets);

            auto it = spurPath.end() - 2;
            if (*it < 0) {
                spurPath.pop_back();
                int toGet = *it * -1;
                spurPath.pop_back();
                for (auto const & vertex : optimalPaths[toGet].second) {
                    spurPath.push_back(vertex);
                }
            }

            delete editableGraph;
            editableGraph = new Graph(graph);
            if (spurPath.empty()) {
                continue;
            }

            auto totalPath = rootPath;
            for (auto node : spurPath) {
                if (node == spurNode) {
                    continue;
                }
                totalPath.push_back(node);
            }
            bool present = false;
            for (auto const & path : candidatePaths) {
                if (totalPath == path.second) {
                    present = true;
                    break;
                }
            }

            if (!present) {
                std::pair<double, std::vector<int>> totalPathPair;
                totalPathPair.first = editableGraph->getPathLength(totalPath, source);
                totalPathPair.second = totalPath;
                candidatePaths.push_back(totalPathPair);
            }
        }

        if (candidatePaths.empty()) {
            break;
        }

        std::sort(candidatePaths.begin(), candidatePaths.end(), PathComparator());

        validPaths.push_back(candidatePaths.front());
        candidatePaths.erase(candidatePaths.begin());

    }
    return validPaths;
}

std::vector<std::pair<double, std::vector<int>>> yen(Graph *graph, int source, int sink, int k) {

    std::vector<std::pair<double, std::vector<int>>> validPaths;
    std::vector<std::pair<double, std::vector<int>>> candidatePaths;

    std::unordered_set<int> visited;

    std::vector<int> optimalPath = dijkstra(graph, source, sink, visited);

    std::pair<double, std::vector<int>> optimalPathPair;
    optimalPathPair.first = graph->getPathLength(optimalPath, source);
    optimalPathPair.second = optimalPath;
    validPaths.push_back(optimalPathPair);

    auto editableGraph = new Graph(graph);
    for (int i=1;i<=k;i++) {
        for (int j=0;j<=validPaths[i-1].second.size() - 2;j++) {
            visited.clear();
            int spurNode = validPaths[i-1].second[j];
            std::vector<int> rootPath = validPaths[i-1].second;
            auto pos = std::find(rootPath.begin(), rootPath.end(), spurNode);
            rootPath.erase(pos+1, rootPath.end());

            for (auto path : validPaths) {
                auto newPos = std::find(path.second.begin(), path.second.end(), spurNode);
                path.second.erase(newPos+1, path.second.end());
                if (rootPath == path.second) {
                    int start = *newPos;
                    int dest = *(newPos+1);
                    editableGraph->removeEdge(start, dest);
                }
            }

            for (auto vertex : rootPath) {
                if (vertex == spurNode) {
                    break;
                }
                visited.insert(vertex);
            }

            auto spurPath = dijkstra(editableGraph, spurNode, sink, visited);
            delete editableGraph;
            editableGraph = new Graph(graph);
            if (spurPath.empty()) {
                continue;
            }

            auto totalPath = rootPath;
            for (auto node : spurPath) {
                if (node == spurNode) {
                    continue;
                }
                totalPath.push_back(node);
            }
            bool present = false;
            for (auto const & path : candidatePaths) {
                if (totalPath == path.second) {
                    present = true;
                    break;
                }
            }

            if (!present) {
                std::pair<double, std::vector<int>> totalPathPair;
                totalPathPair.first = editableGraph->getPathLength(totalPath, source);
                totalPathPair.second = totalPath;
                candidatePaths.push_back(totalPathPair);
            }
        }

        if (candidatePaths.empty()) {
            break;
        }

        std::sort(candidatePaths.begin(), candidatePaths.end(), PathComparator());

        validPaths.push_back(candidatePaths.front());
        candidatePaths.erase(candidatePaths.begin());

    }
    return validPaths;
}

void printPath(std::vector<int> const & path) {
    for (auto node : path) {
        std::cout << "->" << node;
    }
    std::cout << std::endl;
}
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
    auto *reverseGraph = new Graph(vertices, edges);

    for (int i=0;i<edges;i++) {
        int start, end;
        double weight;
        input >> start;
        input >> end;
        input >> weight;
        if (originalGraph->getEdge(start, end) == Edge(-1, -1, 0)) {
            originalGraph->insertEdge(start, end, weight);
        }
        if (reverseGraph->getEdge(end, start) == Edge(-1, -1, 0)) {
            reverseGraph->insertEdge(end, start, weight);
        }
    }

    int source, destination, k;

    input >> source;
    input >> destination;
    input >> k;

    int num = 0;

    /*auto startYen = std::chrono::system_clock::now();
    auto pathPairs = yen(originalGraph, source, destination, k-1);
    auto endYen = std::chrono::system_clock::now();
    std::chrono::duration<double> yenTimeTaken = endYen - startYen;*/

    auto startNewYen = std::chrono::system_clock::now();
    auto result = newYen(originalGraph, reverseGraph, source, destination, k-1);
    auto endNewYen = std::chrono::system_clock::now();
    std::chrono::duration<double> newYenTimeTaken = endNewYen - startNewYen;

    /*std::cout << yenTimeTaken.count() << std::endl;
    for (auto const & pathPair : pathPairs) {
        num++;
        if (num == 0) {
            continue;
        }
        std::cout << "Yen K = " << num << ": " << pathPair.first << " ";
        printPath(pathPair.second);
        std::cout << std::endl;
    }*/

    std::cout << newYenTimeTaken.count() << std::endl;
    num = 0;
    for (auto const & pathPair : result) {
        num++;
        if (num == 0) {
            continue;
        }
        std::cout << "New Yen K = " << num << ": " << pathPair.first << " ";
        printPath(pathPair.second);
        std::cout << std::endl;
    }

    delete originalGraph;
    delete reverseGraph;

    return 0;
}