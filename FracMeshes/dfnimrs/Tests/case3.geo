// Gmsh project created on Tue Jun 14 15:19:12 2022
SetFactory("OpenCASCADE");
gridsize =1;
nx=4;
ny=9;
nz=6;
nstrud = 2;

Point(1)={0.0,0.0,0.0,gridsize};
Point(2)={1,0.0,0.0,gridsize};
Point(3)={1,2.25,0.0,gridsize};
Point(4)={0,2.25,0.0,gridsize};

Line(5)={1,2};
Line(6)={2,3};
Line(7)={3,4};
Line(8)={4,1};

Line Loop(9) ={5,6,7,8};
Plane Surface(10) =9;

Transfinite Line{5,6,7,8} = ny+1;
Transfinite Line{5,7} = nx+1;
Transfinite Surface{10};
Recombine Surface{10};

newEntitiess2[0]=10;
For w In {0:nstrud}
	newEntitiess []=
	Extrude{0,0,1.0/3.0}
	{
		Surface{newEntitiess2[0]};
		Layers{nz/(nstrud+1)};
		Recombine;
	};
	newEntitiess2 = newEntitiess;
EndFor
	
Physical Surface("noflux") = {10, 25, 18, 21, 11, 14, 19, 24, 22, 17, 12};
Physical Surface("inlet") = {16};
Physical Surface("outlet") = {13,23};
Physical Volume("k33") = {3, 2, 1};
