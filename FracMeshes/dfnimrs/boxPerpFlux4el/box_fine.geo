//  Geo file generated by DFNMesh project
//	https://github.com/labmec/dfnMesh// POINTS DEFINITION 

h = 0;

Point(1) = {-1,-1,1, h};
Point(2) = {-1,0,1, h};
Point(3) = {-1,0,-1, h};
Point(4) = {-1,-1,-1, h};
Point(5) = {0,0,1, h};
Point(6) = {0,-1,1, h};
Point(7) = {0,0,-1, h};
Point(8) = {0,-1,-1, h};
Point(9) = {-1,1,1, h};
Point(10) = {-1,1,-1, h};
Point(11) = {0,1,1, h};
Point(12) = {0,1,-1, h};
Point(13) = {1,-1,-1, h};
Point(14) = {1,-1,1, h};
Point(15) = {1,0,-1, h};
Point(16) = {1,0,1, h};
Point(17) = {1,1,-1, h};
Point(18) = {1,1,1, h};


// LINES DEFINITION 

Line(21) = {1,4};
Line(22) = {4,3};
Line(23) = {2,3};
Line(24) = {2,1};
Line(25) = {6,5};
Line(26) = {1,6};
Line(27) = {2,5};
Line(28) = {4,8};
Line(29) = {8,7};
Line(30) = {3,7};
Line(31) = {8,6};
Line(32) = {3,10};
Line(33) = {10,9};
Line(34) = {9,2};
Line(35) = {9,11};
Line(36) = {11,5};
Line(37) = {10,12};
Line(38) = {12,11};
Line(39) = {12,7};
Line(40) = {6,14};
Line(41) = {14,13};
Line(42) = {13,8};
Line(43) = {13,15};
Line(44) = {7,15};
Line(45) = {14,16};
Line(46) = {15,16};
Line(47) = {5,16};
Line(48) = {15,17};
Line(49) = {17,12};
Line(50) = {17,18};
Line(51) = {18,11};
Line(52) = {16,18};
Line(53) = {7,5};


// FACES DEFINITION 

Curve Loop(1) = {-24,23,-22,-21};
Surface(1) = {1};
Curve Loop(2) = {-27,24,26,25};
Surface(2) = {2};
Curve Loop(3) = {30,-29,-28,22};
Surface(3) = {3};
Curve Loop(4) = {28,31,-26,21};
Surface(4) = {4};
Curve Loop(5) = {-34,-33,-32,-23};
Surface(5) = {5};
Curve Loop(6) = {34,27,-36,-35};
Surface(6) = {6};
Curve Loop(7) = {33,35,-38,-37};
Surface(7) = {7};
Curve Loop(8) = {37,39,-30,32};
Surface(8) = {8};
Curve Loop(9) = {-42,-41,-40,-31};
Surface(9) = {9};
Curve Loop(10) = {44,-43,42,29};
Surface(10) = {10};
Curve Loop(11) = {41,43,46,-45};
Surface(11) = {11};
Curve Loop(12) = {-47,-25,40,45};
Surface(12) = {12};
Curve Loop(13) = {-49,-48,-44,-39};
Surface(13) = {13};
Curve Loop(14) = {38,-51,-50,49};
Surface(14) = {14};
Curve Loop(15) = {-46,48,50,-52};
Surface(15) = {15};
Curve Loop(16) = {36,47,52,51};
Surface(16) = {16};
Curve Loop(54) = {25,-53,-29,31};
Surface(54) = {54};
Curve Loop(55) = {30,53,-27,23};
Surface(55) = {55};
Curve Loop(56) = {53,-36,-38,39};
Surface(56) = {56};
Curve Loop(57) = {46,-47,-53,44};
Surface(57) = {57};



// VOLUMES DEFINITION 

Surface Loop(2) = {1,2,3,4,54,55};
Volume(2) = {2};
Surface Loop(3) = {5,6,7,8,55,56};
Volume(3) = {3};
Surface Loop(4) = {9,10,11,12,54,57};
Volume(4) = {4};
Surface Loop(5) = {13,14,15,16,56,57};
Volume(5) = {5};


// COARSE ELEMENTS GROUPING

Physical Volume("c0",17) = {2};
Physical Volume("c1",18) = {3};
Physical Volume("c2",19) = {4};
Physical Volume("c3",20) = {5};



 // FRACTURES

frac0[] = {54,56};

Physical Surface("Fracture0",300) = {frac0[]};



// BOUNDARY CONDITIONS

Physical Surface("inlet",2) = {1,5};
Physical Surface("outlet",3) = {11,15};
Physical Surface("noflux",4) = {2,3,4,6,7,8,9,10,12,13,14,16};


BCfrac0[] = { 25,29,31,36,38,39};
Physical Curve("BCfrac0", 301) = {BCfrac0[]};

// INTER-FRACTURE INTERSECTIONS


// OPTIONS

Coherence Mesh;
Transfinite Curve {:} = 2;
Transfinite Surface{:};
Transfinite Volume{:};
Recombine Surface{:};
Recombine Volume{:};