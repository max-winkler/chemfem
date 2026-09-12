// Channel with a circular obstacle, the geometry of the DFG flow around a cylinder
// benchmark: 2.2 x 0.41, cylinder of radius 0.05 centered at (0.2, 0.2).
//
//   gmsh -2 meshes/cylinder.geo -o meshes/cylinder.msh

L = 2.2;    H = 0.41;
cx = 0.2;   cy = 0.2;   r = 0.05;

h_channel = 0.04;    // element size at the walls
h_cylinder = 0.012;  // element size at the cylinder

Point(1) = {0, 0, 0, h_channel};
Point(2) = {L, 0, 0, h_channel};
Point(3) = {L, H, 0, h_channel};
Point(4) = {0, H, 0, h_channel};

Line(1) = {1, 2};   // bottom wall
Line(2) = {2, 3};   // outflow
Line(3) = {3, 4};   // top wall
Line(4) = {4, 1};   // inflow

Point(5) = {cx, cy, 0, h_cylinder};
Point(6) = {cx-r, cy, 0, h_cylinder};
Point(7) = {cx, cy-r, 0, h_cylinder};
Point(8) = {cx+r, cy, 0, h_cylinder};
Point(9) = {cx, cy+r, 0, h_cylinder};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};
