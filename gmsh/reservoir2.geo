// Gmsh project created on Sun Feb 23 16:19:40 2020
SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1.0, 0, 0, 1.0};
Point(3) = {1.0, 1.0, 0, 1.0};
Point(4) = {0.0, 1.0, 0, 1.0};
Point(5) = {0.5, 0.5, 0, 1.0};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Point(6) = {0.7, 0.5, 0, 1.0};
Point(7) = {0.5, 0.7, 0, 1.0};
Point(8) = {0.3, 0.5, 0, 1.0};
Point(9) = {0.5, 0.3, 0, 1.0};



Point(10) = {0, 0, 1, 1.0};
Point(11) = {1.0, 0, 1, 1.0};
Point(12) = {1.0, 1.0, 1, 1.0};
Point(13) = {0.0, 1.0, 1, 1.0};
Point(14) = {0.5, 0.5, 1, 1.0};

Line(5)={10,11};
Line(6)={11,12};
Line(7)={12,13};
Line(8)={13,10};

Point(15) = {0.7, 0.5, 1, 1.0};
Point(16) = {0.5, 0.7, 1, 1.0};
Point(17) = {0.3, 0.5, 1, 1.0};
Point(18) = {0.5, 0.3, 1, 1.0};//+
Line(9) = {10, 1};
//+
Line(10) = {11, 2};
//+
Line(11) = {12, 3};
//+
Line(12) = {13, 4};
//+
Curve Loop(1) = {4, -9, -8, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 10, -1, -9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, 11, -2, -10};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, 12, -3, -11};
//+
Plane Surface(4) = {4};
//+
Circle(13) = {17, 14, 18};
//+
Circle(14) = {18, 14, 15};
//+
Circle(15) = {15, 14, 16};
//+
Circle(16) = {16, 14, 17};
//+
Curve Loop(5) = {8, 5, 6, 7};
//+
Curve Loop(6) = {15, 16, 13, 14};
//+
Plane Surface(5) = {5, 6};
//+
Circle(17) = {5, 9, 6};
//+
Circle(17) = {6, 5, 9};
//+
Circle(18) = {9, 5, 8};
//+
Circle(19) = {8, 5, 7};
//+
Circle(20) = {7, 5, 6};
//+
Line(21) = {6, 15};
//+
Line(22) = {9, 18};
//+
Line(23) = {8, 17};
//+
Line(24) = {7, 16};
//+
Curve Loop(7) = {19, 24, 16, -23};
//+
Plane Surface(6) = {7};
//+
Curve Loop(8) = {19, 24, 16, -23};
//+
Surface(6) = {8};
//+
Curve Loop(10) = {20, 21, 15, -24};
//+
Surface(7) = {10};
//+
Curve Loop(12) = {17, 22, 14, -21};
//+
Surface(8) = {12};
//+
Curve Loop(14) = {18, 23, 13, -22};
//+
Surface(9) = {14};
//+
Curve Loop(16) = {2, 3, 4, 1};
//+
Curve Loop(17) = {18, 19, 20, 17};
//+
Plane Surface(10) = {16, 17};
//+
Surface Loop(1) = {10, 3, 5, 1, 2, 4, 7, 6, 9, 8};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {7, 6, 9, 5, 1, 2, 3, 4, 8, 10};
//+
Physical Surface("Well") = {6, 9, 7, 8};
//+
Physical Surface("Boundary") = {4, 5, 2, 3, 10, 1};
//+
Physical Volume("Volumen") = {1};
