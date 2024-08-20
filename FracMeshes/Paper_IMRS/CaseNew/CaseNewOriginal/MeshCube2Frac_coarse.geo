// Gmsh project created on Thu Jan 18 19:10:32 2024
SetFactory("OpenCASCADE");
h=1;
xmax=2;
ymax=2;
zmax=1;

Point(1)={0,0,0,h};
Point(2)={0,ymax,0,h};
Point(3)={xmax,ymax,0,h};
Point(4)={xmax,0,0,h};

Point(5)={0,0,zmax,h};
Point(6)={0,ymax,zmax,h};
Point(7)={xmax,ymax,zmax,h};
Point(8)={xmax,0,zmax,h};//+


Point(9)={0,ymax/2,zmax,h};
Point(10)={xmax,ymax/2,zmax,h};

Point(11)={0,ymax/2,0,h};
Point(12)={xmax,ymax/2,0,h};



Point(13)={xmax/2,ymax,0,h};
Point(14)={xmax/2,0,0,h};

Point(15)={xmax/2,ymax,zmax,h};
Point(16)={xmax/2,0,zmax,h};

Point(17)={xmax/2,ymax/2,zmax,h};
Point(18)={xmax/2,ymax/2,0,h};



//+
Line(1) = {5, 16};
//+
Line(2) = {16, 8};
//+
Line(3) = {8, 4};
//+
Line(4) = {4, 14};
//+
Line(5) = {14, 1};
//+
Line(6) = {1, 5};
//+
Line(7) = {9, 17};
//+
Line(8) = {17, 10};
//+
Line(9) = {10, 12};
//+
Line(10) = {12, 18};
//+
Line(11) = {18, 11};
//+
Line(12) = {11, 9};
//+
Line(13) = {7, 3};
//+
Line(14) = {3, 13};
//+
Line(15) = {13, 2};
//+
Line(16) = {2, 6};
//+
Line(17) = {6, 15};
//+
Line(18) = {15, 7};
//+
Line(19) = {16, 14};
//+
Line(20) = {17, 18};
//+
Line(21) = {15, 13};
//+
Line(22) = {5, 9};
//+
Line(23) = {1, 11};
//+
Line(24) = {11, 2};
//+
Line(25) = {9, 6};
//+
Line(26) = {18, 13};
//+
Line(27) = {17, 15};
//+
Line(28) = {14, 18};
//+
Line(29) = {16, 17};
//+
Line(30) = {12, 3};
//+
Line(31) = {4, 12};
//+
Line(32) = {10, 7};
//+
Line(33) = {8, 10};
//+
Curve Loop(1) = {1, 19, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 4, -19};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, 20, 11, 12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, 9, 10, -20};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {17, 21, 15, 16};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {18, 13, 14, -21};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {25, -16, -24, 12};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {22, -12, -23, 6};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {29, 20, -28, -19};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {27, 21, -26, -20};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {32, 13, -30, -9};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {9, -31, -3, 33};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {5, 23, -11, -28};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {1, 29, -7, -22};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {2, 33, -8, -29};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {4, 28, -10, -31};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {11, 24, -15, -26};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {7, 27, -17, -25};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {10, 26, -14, -30};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {8, 32, -18, -27};
//+
Plane Surface(20) = {20};
//+
Surface Loop(1) = {14, 1, 13, 8, 3, 9};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {16, 2, 15, 12, 9, 4};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {19, 6, 20, 11, 10, 4};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {5, 18, 7, 17, 10, 3};
//+
Volume(4) = {4};
//+
Physical Surface("inlet", 34) = {7, 8, 1};
//+
Physical Surface("outlet", 35) = {11};
//+
Physical Surface("noflux", 36) = {5, 6, 19, 20, 18, 17, 14, 13, 2, 16, 15, 12};
//+
Physical Volume("k11", 37) = {4, 3, 1, 2};

Coherence Mesh;
Transfinite Curve {:} = 2;
Transfinite Surface{:};
Transfinite Volume{:};
Recombine Surface{:};
Recombine Volume{:};
