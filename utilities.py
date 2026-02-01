import random as r 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.colors import ListedColormap

class Atom: 
    def __init__(self, id: int, x: int, y: int):
        self.id = id 
        self.x = x
        self.y = y
        self.fixed = False

    def x_shift(self, threshold: float):
        """
        Moves the atom left or right.
        If a random roll is greater than the threshold, move in a random direction.
        """
        # Generate a value between 0 and 1 to check against threshold
        if r.random() < threshold:
            # Choose move: -1 (left) or 1 (right)
            move = r.choice([-1, 1])
            self.x += move

    def y_shift(self, threshold: float):
        """
        Moves the atom up or down.
        If a random roll is greater than the threshold, move in a random direction.
        """
        # Generate a value between 0 and 1 to check against threshold
        if r.random() < threshold:
            # Choose move: -1 (left) or 1 (right)
            move = r.choice([-1, 1])
            self.y += move

    def xy_shift(self, threshold: float): 
        direction = r.choice(['x', 'y'])
        if direction =='x': 
            self.x_shift(threshold)
        elif direction == 'y': 
            self.y_shift(threshold)

    def get_atom_pos(self): 
        return((self.x, self.y))

    def set_atom_pos(self, pos:tuple):
        self.x = pos[0]
        self.y = pos[1]

    def __repr__(self):
        return f"Atom(id={self.id}, x={self.x}, y={self.y}, fixed={self.fixed})"

class Layer: 
    def __init__(self, id:int, sites_1d: int):
        self.id = id
        self.atoms = []
        self.sites_1d = sites_1d
        self.max_sites = self.sites_1d**2
        self.n_atoms = len(self.atoms)
        self.atom_map = {}

    @property
    def occ_sites(self):
        """Automatically returns a fresh list of (x, y) tuples from current atoms."""
        return [atom.get_atom_pos() for atom in self.atoms]

    def __repr__(self):
        return f"Layer(id={self.id}, n_atoms={len(self.atoms)})"
        pass

    def add_atom(self, atom: Atom): 
        pos = atom.get_atom_pos()
        if pos not in self.atom_map:
            self.atoms.append(atom) 
            self.atom_map[pos] = atom

    def get_neighbors(self, pos):
        """Helper to find atoms in adjacent cardinal directions."""
        count = 0
        directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        for dx, dy in directions:
            if (pos[0] + dx, pos[1] + dy) in self.atom_map:
                count += 1
        return count
    
    def add_rand_atom(self):
        added_x = r.randint(0, self.sites_1d-1)
        added_y = r.randint(0, self.sites_1d-1)
        pos = (added_x, added_y)
        if not (added_x, added_y) in self.occ_sites:
            added_atom = Atom(len(self.occ_sites), added_x, added_y)
            self.atoms.append(added_atom) 
            self.occ_sites.append(added_atom.get_atom_pos())
            self.atom_map[pos] = added_atom

    def rand_walk(self):
        # We iterate over a copy of the list to ensure stability
        for atom in self.atoms:
            if atom.fixed:
                continue

            # 1. Update neighbor count and calculate move probability
            # Probability formula: 1.0 (no neighbors) down to 0.0 (4 neighbors)
            atom.neighbor_count = self.get_neighbors(atom.get_atom_pos())
            
            # Linear decrease: Each neighbor reduces move chance by 25%
            move_probability = max(0.0, 1.0 - (0.6+ atom.neighbor_count * 0.1))

            if move_probability <= 0:
                atom.fixed = True
                continue

            # 2. Attempt Move
            old_pos = atom.get_atom_pos()
            atom.xy_shift(move_probability)
            new_pos = atom.get_atom_pos()

            # 3. Validation Logic
            out_of_bounds = (
                new_pos[0] < 0 or new_pos[0] >= self.sites_1d or
                new_pos[1] < 0 or new_pos[1] >= self.sites_1d
            )
            collision = new_pos in self.atom_map and self.atom_map[new_pos] != atom

            if out_of_bounds or collision:
                atom.set_atom_pos(old_pos)
            else:
                # Update map if position actually changed
                if old_pos != new_pos:
                    del self.atom_map[old_pos]
                    self.atom_map[new_pos] = atom
            
            # 4. Final state update: If surrounded after moving, fix it
            atom.neighbor_count = self.get_neighbors(atom.get_atom_pos())
            if atom.neighbor_count >= 4:
                atom.fixed = True


    def plot_layer(self):

        grid = np.zeros((self.sites_1d, self.sites_1d))
        for x, y in self.occ_sites: 
            grid[y, x] = 1 

        cmap = ListedColormap(['white', 'blue'])

        plt.imshow(grid, cmap=cmap, origin='lower', extent=[0, self.sites_1d, 0, self.sites_1d])

        # Optional: Add grid lines and labels for better visibility
        plt.grid(True, which='both', color='lightgrey', linestyle='-', linewidth=1)
        plt.xticks(np.arange(0, self.sites_1d + 1, 1))
        plt.yticks(np.arange(0, self.sites_1d + 1, 1))
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title(f"{self.sites_1d}x{self.sites_1d} Grid: Occupied Sites (Blue)")

        plt.show()


    








