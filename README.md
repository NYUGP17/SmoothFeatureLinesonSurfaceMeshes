# Smooth Feature Lines on Surface Meshes

## Goal of this project

In this project, the author will implement an algorithm to extract feature lines of 3D objects. The algorithm is described in [Smooth Feature Lines on Surface Meshes](http://ddg.math.uni-goettingen.de/pub/feature.pdf). The author will construct a discrete differential geometry and use a filter to improve the stability and smoothness of the extracted lines.

## Tasks

### Discrete extremality

- Calculate the mean curvature and the shape operator at each edge e
- Calculate the vertex-based shape operator and the discrete principal curvatures
- Make consistent choice of sign of K_i
- Remove singular meshes


### Smooth discrete extremalities

#### A modification of Laplacian smoothing

- Compute the extremalities at each vertex of the mesh by using an arbitrary choice of sign of K_i
- Compute the modified Laplacian for e_i using the cotan weights.

#### Spatial Fairing

A modification of the smoothing scheme presented in [Anisotropic Filtering of Non-Linear Surface Features](https://graphics.tudelft.nl/~klaus/papers/smooth.pdf). (What is that?)

### Trace feature lines in regular triangles

Compute ridge line segment for each regular triangle

### Process singular triangles

Consider the adjacect regular triangles. 
Mark the endpoint of feature line segment in the common edge.

- Connect the endpoints, if there are 2 marked edges.
- Add a point: barycenter, and connect it with all endpoints, if there are 3 marked edges.
- Do nothing otherwise.

### Remove small ridges by a threshold filter

Discard small and squiggly lines

## Input

A 3D object in certain format (.off).

## Output

Visualization of feature lines.
