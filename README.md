# Different-solutions-to-TSP
This is a coursework project on creating different solutions to the Traveling Salesperson Problem:

NOTE: Copying this code is against University of Michigan's Honor Code. For your convenience here is the link to the honor code: http://ossa.engin.umich.edu/wp-content/uploads/sites/212/2015/04/Honor-code-pamphlet-Adobe-Prof.pdf. 

3 parts to the project:

1. Part A (--mode MST): Finding and printing the Minimum-Spanning-Tree given a graph in the form of vertices (coordinates from an input file)
2. Part B (--mode FASTTSP): Finding a fast and good (not the best) solution to the Shortest Hamiltonian Cycle problem (Traveling Salesperson Problem) using a TSP heuristic called arbitrary insertion/random insertion
  - I selected this heuristic out of all greedy ones possible as it was by-far the fastest and cleanest to implement
  - Previous attempts include: - Nearest Neighbour heuristic with 2-OPT
                               - Cheapest Insertion
                               - Farthest Insertion
3. Part C (--mode OPTTSP): Finding the optimal solution to the Shortest Hamiltonian Cycle problem (Traveling Salesperson Problem):
Here, I use a Branch-and-Bound approach where:
Step 1: Upper-bound (best guess estimate) is provided by running arbitrary insertion
Step 2: genPerms() function generates all permutations but prunes the ones which are NOT promising (decided by the promising function)
        - Promising() : finds the projected potential cost of a certain branch by finding an MST for the remaining vertices in the path that are NOT a part of the partial soln. already

File breakdown:

- GameADT.h holds all code for 3 different parts of the project
- amongus.cpp hold main and handles command line arguments using getopt_long()
- sample-ab/c/d/e/f.txt are sample inputs
- test-1 to test-10 are test files I created for bug-catching

This project was part 4 of the Data Structures and Algorithms course at the University of Michigan. 

NOTE: Among Us is a theme given to the project. It is only used in Part A where zones are assigned to different areas of the cartesian graph
                                                - Outer zone: QI QII and QIV
                                                - Decontamination zone: -ve x-axis and -ve y-axis points
                                                - Lab zone: QIII
      The idea is that one cannot travel from a point in the Outer zone to the Inner zone without traversing through the Decontamination zone
