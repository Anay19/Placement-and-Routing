# Placement-and-Routing

Cell Placement Algorithm

The cell placement algorithm used is the Force-Directed Placement Algorithm. This algorithm aims to place the cells in a way that minimizes wire length and congestion while considering the connectivity between the cells.

Data Structures
    Net Information Vector: Stores information about each net, including the net number, source terminal, target terminal, and source and target cell numbers.
    Connectivity Matrix: Stores the connections between cells through nets.
    Cell Weight Vector: Stores the connectivity size (number of connections) of each cell.
    Cell Placement Vector: Stores the location of each cell in the placement matrix.
    Placement Matrix: Used to manipulate the cells and place them at the best possible locations using the force-directed algorithm.

Algorithm Steps
    Initially, all cells are placed randomly in the placement matrix.
    Cells are sorted based on their connectivity (number of connections).
    The cell with the highest connectivity is chosen as the seed cell.
    The target location for the seed cell is calculated using the zero-force equation, which equates the net force on the cell to zero and finds the corresponding new x and y coordinates.
    The algorithm then proceeds based on the target location:
        If the target location is empty, the seed cell is placed there, and the location is marked as locked.
        If the target location is occupied, the seed cell is moved there, and the occupying cell becomes the new seed.
        If the target location is blocked, the algorithm finds the nearest empty cell and proceeds accordingly.
        If no suitable location is found after a certain number of attempts (abort limit), all locked cells are unlocked, and the process is repeated.
    The algorithm continues until the iteration limit is reached.

The time complexity of the force-directed placement algorithm depends on the number of cells and the abort limit, but it is generally efficient for moderate-sized designs.
Routing Algorithm (Lee Algorithm)

The routing algorithm used is the Lee Algorithm, which is a wave propagation algorithm that finds the shortest path between a source and a target while considering blockages and vias (layer changes).

Data Structures
    Struct for Net Information: Stores information about Lee numbers, blocks, vias, net numbers, and net surroundings in both layers.
    Cell Coordinate Vector: Stores the locations of cell terminals.
    Queue: Used for breadth-first search during wave propagation.

Algorithm Steps
    Wave Propagation: The Lee algorithm propagates a wave from the source to the target by assigning a number equivalent to the Manhattan distance at each block. The wave propagation ends when the target point is reached. This process happens in two layers, and if there is a blockage in one layer, the algorithm checks the other layer and assigns a Lee number if there is no blockage.
    Back Trace: The back trace starts at the target and tries to follow the same direction as the wave propagation until it cannot proceed further. When it encounters a turn, it takes the path with the smallest Lee number compared to its previous value. Whenever there is a layer change (from layer 1 to layer 2 or vice versa), a via is inserted at that location.

The time complexity of the Lee Algorithm is O(L^2), where L is the Manhattan distance between the source and target, for wave propagation, and O(L) for back tracing.

Both the cell placement and routing algorithms use efficient data structures and algorithms to optimize the chip design process. The force-directed placement algorithm aims to minimize wire length and congestion, while the Lee Algorithm finds the shortest path between cells while considering blockages and layer changes.
