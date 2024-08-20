// Gmsh project created on Thu Jan 18 19:10:32 2024
SetFactory("OpenCASCADE");
h=1;
xmax=1;
ymax=1;
zmax=0.25;

Point(1)={0,0,0,h};
Point(2)={0,ymax,0,h};
Point(3)={xmax,ymax,0,h};
Point(4)={xmax,0,0,h};

Point(5)={0,0,zmax,h};
Point(6)={0,ymax,zmax,h};
Point(7)={xmax,ymax,zmax,h};
Point(8)={xmax,0,zmax,h};//+
Line(1) = {2, 6};
//+
Line(2) = {6, 7};
//+
Line(3) = {7, 3};
//+
Line(4) = {3, 2};
//+
Line(5) = {1, 5};
//+
Line(6) = {5, 8};
//+
Line(7) = {8, 4};
//+
Line(8) = {4, 1};
//+
Line(9) = {4, 3};
//+
Line(10) = {8, 7};
//+
Line(11) = {1, 2};
//+
Line(12) = {5, 6};
//+
Curve Loop(1) = {11, 1, -12, -5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -9, -7, 10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, -4, -9, 8};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {2, -10, -6, 12};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 1, 2, 3};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {8, 5, 6, 7};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {4, 5, 3, 1, 6, 2};
//+
Volume(1) = {1};
//+
Physical Surface("bcInlet", 13) = {1};
//+
Physical Surface("bcOutlet", 14) = {2};
//+
Physical Surface("bcNoFlux", 15) = {4, 3, 5, 6};
//+
Physical Volume("k11", 16) = {1};

Coherence Mesh;
Transfinite Curve {1,5,3,7} = 6;
Transfinite Curve {2,4,6,8} = 21;
Transfinite Curve {12,11,9,10} = 6;
Transfinite Surface{:};
Transfinite Volume{:};
Recombine Surface{:};
Recombine Volume{:};
