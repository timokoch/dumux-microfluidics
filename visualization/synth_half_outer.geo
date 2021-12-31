//+
SetFactory("OpenCASCADE");

cylRadius = 1.5;
y = 8.1;

Cylinder(1) = {12.25, y, 0, 0, 0, 6.75, cylRadius, 2*Pi};
//+
Cylinder(2) = {12.25, -y, 0, 0, 0, 6.75, cylRadius, 2*Pi};

h = 5.2;
c = 2*y;
R = 0.5*h + c*c/(8*h);
cl__ = 0.5;
x = R - h;
norm = Sqrt(x*x + y*y);

// circle mid-points for curve
Point(105) = {12.25 - R + h, 0.0, 0.0, cl__};

// thickness
thickness = 0.95*cylRadius;

// inner
Point(106) = {12.25 - thickness*x/norm, y - thickness*y/norm, 0.0, cl__};
Point(107) = {12.25 - thickness*x/norm, -y + thickness*y/norm, 0.0, cl__};
// outer
Point(108) = {12.25 + thickness*x/norm, y + thickness*y/norm, 0.0, cl__};
Point(109) = {12.25 + thickness*x/norm, -y - thickness*y/norm, 0.0, cl__};
Circle(13) = {106, 105, 107};
Circle(14) = {108, 105, 109};
// sides
Line(15) = {106, 108};
Line(16) = {107, 109};

Curve Loop(13) = {13, 16, -14, -15};
Plane Surface(13) = {13};

Extrude {0, 0, 4.75} { Surface{13}; }
BooleanUnion{ Volume{3}; Delete; }{ Volume{1}; Volume{2}; Delete; }

Point(200) = {0.5*(3.462435334542 + 3.798170069313) + 7.3, 0.5*(-6.14367530233 -7.067145485502) - 2, 0.0, cl__};
Point(201) = {0.5*(3.462435334542 + 3.798170069313) + 7.3, -(0.5*(-6.14367530233 -7.067145485502) - 2), 0.0, cl__};
