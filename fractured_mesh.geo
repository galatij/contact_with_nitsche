Mesh.MshFileVersion = 2.2;
Mesh.Format = 1;

// ------------ Parameters ------------
Lx = 1;
Ly = 4;
Lz = 4;
faultX0 = 0.5;
faultX1 = 0.639683;
h = 0.25;

// ------------ Points ------------
Point(1) = {0, 0, 0, h};
Point(2) = {faultX0, 0, 0, h};
Point(3) = {Lx, 0, 0, h};
Point(4) = {0, Ly, 0, h};
Point(5) = {faultX0, Ly, 0, h};
Point(6) = {Lx, Ly, 0, h};
Point(7) = {0, 0, Lz, h};
Point(8) = {faultX1, 0, Lz, h};
Point(9) = {Lx, 0, Lz, h};
Point(10) = {0, Ly, Lz, h};
Point(11) = {faultX1, Ly, Lz, h};
Point(12) = {Lx, Ly, Lz, h};

// ------------ Fault Surface ------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 6};
Line(4) = {6, 5};
Line(5) = {5, 4};
Line(6) = {4,1};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,12};
Line(10) = {12,11};
Line(11) = {11,10};
Line(12) = {10,7};
Line(13) = {1,7};
Line(14) = {3,9};
Line(15) = {6,12};
Line(16) = {4,10};
Line(17) = {2,5};
Line(18) = {5,11};
Line(19) = {11,8};
Line(20) = {8,2};

// ------------ Left Block Faces ------------
Line Loop(1) = {-6,-5,-17,-1};
Plane Surface(1) = {1};
Line Loop(2) = {5,16,-11,-18};
Plane Surface(2) = {2};
Line Loop(3) = {7,-19,11,12};
Plane Surface(3) = {3};
Line Loop(4) = {1,-20,-7,-13};
Plane Surface(4) = {4};
Line Loop(9) = {6,13,-12,-16};
Plane Surface(9) = {9};

// ------------ Right Block Faces ------------
Line Loop(5) = {-2,17,-4,-3};
Plane Surface(5) = {5};
Line Loop(6) = {4,18,-10,-15};
Plane Surface(6) = {6};
Line Loop(7) = {8,9,10,19};
Plane Surface(7) = {7};
Line Loop(8) = {2,14,-8,20};
Plane Surface(8) = {8};
Line Loop(10) = {3,15,-9,-14};
Plane Surface(10) = {10};

// ------------ Fault Surface ------------
Line Loop(11) = {17,18,19,20};
Plane Surface(1000) = {11};

// ------------ Volume Definitions ------------
Surface Loop(1001) = {1,2,3,4,9,1000};
Volume(101) = {1001};
Surface Loop(1002) = {5,6,7,8,10,-1000};
Volume(102) = {1002};

// ------------ Physical Regions ------------
Physical Surface("BottomLeft") = {1};
Physical Surface("YmaxLeft") = {2};
Physical Surface("TopLeft") = {3};
Physical Surface("YminLeft") = {4};
Physical Surface("BottomRight") = {5};
Physical Surface("YmaxRight") = {6};
Physical Surface("TopRight") = {7};
Physical Surface("YminRight") = {8};
Physical Surface("Xmin") = {9};
Physical Surface("Xmax") = {10};
Physical Volume("BulkLeft") = {101};
Physical Volume("BulkRight") = {102};
Physical Surface("Fault") = {1000};

