Ro = 1;
Ri = 0.3;

cx = 0.1;
cy = 0.1;

Point(1) = {0, 0, 0};
Point(2) = {Ro, 0, 0};
Point(3) = {0, Ro, 0};
Point(4) = {-Ro, 0, 0};
Point(5) = {0, -Ro, 0};
Point(6) = {cx, cy, 0};
Point(7) = {cx, cy - Ri, 0};
Point(8) = {cx + Ri, cy, 0};
Point(9) = {cx, cy + Ri, 0};
Point(10) = {cx - Ri, cy, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Circle(5) = {7, 6, 8};
Circle(6) = {8, 6, 9};
Circle(7) = {9, 6, 10};
Circle(8) = {10, 6, 7};

Line Loop(9) = {1, 2, 3, 4, 5, 6, 7 ,8};
Plane Surface(10) = {9};

Transfinite Line{1, 2, 3, 4} = 30;
Transfinite Line{5, 6, 7, 8} = 10;
Recombine Surface{10};

newEntities[] =
Extrude {0,0,0.1}
{
	Surface{10};
	Layers{1};
	Recombine;
};

Physical Surface("backAndFront") = {10,newEntities[0]};
Physical Surface("outerWall") = {newEntities[2],newEntities[3],newEntities[4],newEntities[5]};
Physical Surface("deformedWall") = {newEntities[6], newEntities[7], newEntities[8], newEntities[9]};
Physical Volume(60) = {newEntities[1]};

Mesh 3;
Save "plateWithHole.msh";
