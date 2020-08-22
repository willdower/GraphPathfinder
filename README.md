# GraphPathfinder

This program, given a weighted directed graph, searches it to find the most efficient path from the start to the goal using Dijkstra's algorithm. It then innovates on Yen's to find the k-shortest paths in quicker time.

The innovation works by flipping the direction of every edge and running Dijkstra's backwards. This means that the shortest path from each node to the goal will be found and stored, so that much of the calculation required by Yen's has already been done and is not repeated multiple times, trading short-term speed for long-term efficiency, which is faster on larger graphs.
