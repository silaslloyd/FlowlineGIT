//+
L={L};
GL={GL};
W=1000;
resb={BoundaryRes};
resbb={BoundaryRes};
resgl={GLRes};
Point(1) = {0, 0, 0, resb};
Point(2) = {GL - GL/16, 0, 0, resgl};
Point(3) = {GL, 0, 0, resgl};
Point(4) = {GL+GL/16, 0, 0, resgl};
Point(5) = {L, 0, 0, resb};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
Line(3) = {3,4};
Line(4) = {4,5};
//+
Extrude {0, W, 0} {
  Curve{1}; Curve{2};Curve{3};Curve{4};Recombine;Layers{1};
}
//+


//+
Physical Surface(21) = {8, 12, 16, 20};
//+
Physical Curve(22) = {6};
//+
Physical Curve(23) = {19};
//+
Physical Curve(24) = {1, 2, 3, 4};
//+
Physical Curve(25) = {5, 9, 13, 17};
