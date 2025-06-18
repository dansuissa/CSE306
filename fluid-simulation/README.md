# Labs 6 & 7: Voronoi Diagrams and Optimal Transport

---

## Lab 6: Voronoi Diagrams

This program computes the Voronoi diagram of a set of 2D points using an optimized O(N log N) algorithm accelerated with a k-d tree (`nanoflann`).

### How to Build and Run

```bash
# Build the executable
make voronoi

# Run the program
./voronoi
```

This will generate `cvt_uniform_optimized.svg`, which contains the final Voronoi diagram.

---

## Lab 7: Optimal Transport

This program computes a semi-discrete optimal transport map between a set of weighted Dirac masses and a uniform target density on the unit square. It uses the L-BFGS optimization library to find the optimal weights for the corresponding power diagram.

### How to Build and Run

```bash
# Build the executable
make optimal_transport

# Run the program
./optimal_transport
```

This will generate `power_uniform.svg`, showing a power diagram where cell areas correspond to a centered Gaussian distribution.

---

### Build All

To build both executables at once, simply run:
```bash
make
```

### Clean Up
```bash
make clean
``` 