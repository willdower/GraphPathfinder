#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <limits>
#include <algorithm>

struct path {
    std::vector<int> pathVector;
    double totalWeight;

    path() {
        totalWeight = 0;
    }
};

class pathCompare {
    bool operator () (const path & pathOne, const path & pathTwo) {
        return (pathOne.totalWeight < pathTwo.totalWeight);
    }
};

int main(int argc, char **argv) {

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

    auto **adjacencyMatrix = new double * [vertices];
    for (int i=0;i<vertices;i++) {
        adjacencyMatrix[i] = new double [vertices]{0};
    }

    for (int i=0;i<edges;i++) {
        int start, end;
        double weight;
        input >> start;
        input >> end;
        input >> weight;
        adjacencyMatrix[start][end] = weight;
    }

    int source, destination, k;

    input >> source;
    input >> destination;
    input >> k;

    auto ***paths = new path ** [vertices];
    for (int i=0;i<vertices;i++) {
        paths[i] = new path * [k];
        for (int j=0;j<k;j++) {
            paths[i][j] = new path;
            paths[i][j]->totalWeight = std::numeric_limits<double>::max(); // Set each path's weight to max so it will be lowest priority
        }
    }

    std::vector<int> nodeQueue;
    nodeQueue.push_back(source); // Push first node
    std::vector<int> visitedNodes;
    paths[source][0] = new path;
    paths[source][0]->pathVector.push_back(source);


    while (!nodeQueue.empty()) {

        double currentNodeWeight = paths[nodeQueue[0]][0]->totalWeight, pos = 0;

        for (int i = 0; i < nodeQueue.size(); i++) {
            if (paths[nodeQueue[i]][0]->totalWeight < currentNodeWeight) {
                currentNodeWeight = paths[nodeQueue[i]][0]->totalWeight;
                pos = i;
            }
        }
        int currentNode = nodeQueue[pos];
        nodeQueue.erase(nodeQueue.begin() + pos);

        //std::cout << "Current Node: " << currentNode << std::endl;
        for (int i = 0; i < vertices; i++) { // For each node connected to the current node
            double currentWeight = adjacencyMatrix[currentNode][i]; // Get the connection weight
            if (currentWeight == 0) {
                // No edge
                continue;
            }
            for (int j = 0; j < k; j++) { // For each path to the current node
                //std::cout << "Loop B" << std::endl;
                // Add each path as a path to the new node, with the new connection added
                if (paths[currentNode][j]->totalWeight == std::numeric_limits<double>::max()) {
                    break;
                }
                if (std::find(paths[currentNode][j]->pathVector.begin(), paths[currentNode][j]->pathVector.end(), i) !=
                    paths[currentNode][j]->pathVector.end()) {
                    continue;
                }

                path *newPath = new path;
                newPath->pathVector = paths[currentNode][j]->pathVector;
                newPath->totalWeight = paths[currentNode][j]->totalWeight;
                newPath->pathVector.push_back(i);
                newPath->totalWeight += currentWeight;

                if (i == destination) {
                    std::cout << newPath->totalWeight << std::endl;
                }

                int pos = 0;
                bool insert = true;
                // Check if new path is better than the top-K paths already there
                while (paths[i][pos]->totalWeight < newPath->totalWeight) {
                    //std::cout << "Loop C" << std::endl;
                    pos++;
                    if (pos >= k) {
                        insert = false;
                        break;
                    }
                }
                // If it is to be inserted, push paths behind back and then insert new path
                if (insert) {
                    delete paths[i][k - 1];
                    for (int n = k - 2; n >= pos; n--) {
                        paths[i][n + 1] = paths[i][n];
                    }
                    paths[i][pos] = newPath;
                } else {
                    delete newPath;
                }
            }
            if (std::find(visitedNodes.begin(), visitedNodes.end(), i) ==
                visitedNodes.end()) { // If node hasn't been popped before
                nodeQueue.push_back(i);
                visitedNodes.push_back(i);
            }
        }
    }


    std::cout << "------OUTPUT-----" << std::endl;
    for (int i = 0; i < k; i++) {
        if (paths[destination][i]->totalWeight == std::numeric_limits<double>::max()) {
            break;
        }
        std::cout << "Weight: " << paths[destination][i]->totalWeight << std::endl;
        for (auto j : paths[destination][i]->pathVector) {
            std::cout << "->" << j;
        }
        std::cout << std::endl;
    }

    return 0;
}
