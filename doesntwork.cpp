//
// Created by William on 11/05/2020.
//

#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <limits>
#include <algorithm>
#include <array>

struct node {
    int val;
    int level;
};

struct path {
    std::vector<node*> pathVector;
    double totalWeight = 0;
};

void update(path ***paths, node **existingNodes, double **adjacencyMatrix, int vertices, node *currentNode, int currentLevel, int k, std::queue<node*> & toUpdate) {
    // For each edge connected to the current node
    for (int i=0;i<vertices;i++) {
        if (adjacencyMatrix[currentNode->val][i] == 0) {
            continue;
        }

        if (existingNodes[i] == nullptr) {
            // Allocate destination node if it doesn't exist yet
            existingNodes[i] = new node;
            existingNodes[i]->val = i;
            existingNodes[i]->level = currentLevel;
        }

        for (int j=0;j<k;j++) {
            if (paths[currentNode->val][j]->totalWeight == std::numeric_limits<double>::max()) {
                break;
            }
            path *newPath = new path;
            newPath->pathVector = paths[currentNode->val][j]->pathVector;
            newPath->pathVector.push_back(existingNodes[i]);
            newPath->totalWeight = paths[currentNode->val][j]->totalWeight+adjacencyMatrix[currentNode->val][i];

            int pos = 0;
            bool insert = true;
            // Check if new path is better than the top-K paths already there
            while (paths[i][pos]->totalWeight < newPath->totalWeight && paths[i][pos]->pathVector != newPath->pathVector) {
                pos++;
                if (pos >= k) {
                    // Not better than top k
                    insert = false;
                    break;
                }
            }
            if (paths[i][pos]->pathVector != newPath->pathVector) {
                insert = false;
            }
            // If it is to be inserted, push paths behind back and then insert new path
            if (insert) {
                toUpdate.push(existingNodes[i]);
                delete paths[i][k - 1];
                for (int n = k - 2; n >= pos; n--) {
                    paths[i][n + 1] = paths[i][n];
                }
                paths[i][pos] = newPath;
            } else {
                delete newPath;
            }

        }

    }
}

void updateOnSameLevel(path ***paths, node **existingNodes, double **adjacencyMatrix, std::queue<node*> nodeQueue, int vertices, int currentLevel, int k) {
    int queueSize = nodeQueue.size();
    for (int i=0;i<queueSize;i++) {
        node *currentNode = nodeQueue.front();
        nodeQueue.pop();
        for (int j=0;j < vertices; j++) {
            std::queue<node*> toUpdate;
            update(paths, existingNodes, adjacencyMatrix, vertices, currentNode, currentLevel, k, toUpdate);
            updateOnSameLevel(paths, existingNodes, adjacencyMatrix, toUpdate, vertices, currentLevel, k);
        }
    }
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

    node *sourceNode = new node;
    sourceNode->level = 0;
    sourceNode->val = source;
    paths[sourceNode->val][0]->totalWeight = 0;
    paths[sourceNode->val][0]->pathVector.push_back(sourceNode);

    std::queue<node*> nodeQueue;
    node **existingNodes = new node* [vertices];
    for (int i=0;i<vertices;i++) {
        existingNodes[i] = nullptr;
    }
    nodeQueue.push(sourceNode);
    int currentLevel = 0, membersInLevel = 1;
    std::vector<node*> visitedNodes;

    while (!nodeQueue.empty()) {
        node *currentNode = nodeQueue.front();
        nodeQueue.pop();

        // For each edge connected to the current node
        for (int i=0;i<vertices;i++) {
            if (adjacencyMatrix[currentNode->val][i] == 0) {
                continue;
            }

            if (existingNodes[i] == nullptr) {
                // Allocate destination node if it doesn't exist yet
                existingNodes[i] = new node;
                existingNodes[i]->val = i;
                existingNodes[i]->level = currentLevel;
            }

            for (int j=0;j<k;j++) {
                if (paths[currentNode->val][j]->totalWeight == std::numeric_limits<double>::max()) {
                    break;
                }
                path *newPath = new path;
                newPath->pathVector = paths[currentNode->val][j]->pathVector;
                newPath->pathVector.push_back(existingNodes[i]);
                newPath->totalWeight = paths[currentNode->val][j]->totalWeight+adjacencyMatrix[currentNode->val][i];

                int pos = 0;
                bool insert = true;
                // Check if new path is better than the top-K paths already there
                while (paths[i][pos]->totalWeight < newPath->totalWeight) {
                    pos++;
                    if (pos >= k) {
                        // Not better than top k
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

            if (std::find(visitedNodes.begin(), visitedNodes.end(), existingNodes[i]) == visitedNodes.end()) {
                nodeQueue.push(existingNodes[i]);
                visitedNodes.push_back(existingNodes[i]);
            }
        }

        membersInLevel--;
        if (membersInLevel == 0) {
            membersInLevel = nodeQueue.size();
            currentLevel++;
            updateOnSameLevel(paths, existingNodes, adjacencyMatrix, nodeQueue, vertices, currentLevel, k);
        }

    }

    std::cout << "------OUTPUT-----" << std::endl;
    for (int i = 0; i < k; i++) {
        if (paths[destination][i]->totalWeight == std::numeric_limits<double>::max()) {
            break;
        }
        std::cout << "Weight: " << paths[destination][i]->totalWeight << std::endl;
        for (auto j : paths[destination][i]->pathVector) {
            std::cout << "->" << j->val;
        }
        std::cout << std::endl;
    }

    return 0;
}