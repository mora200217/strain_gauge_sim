SetFactory("OpenCASCADE");
Merge "../step/load_cell_beam.step";
//+
Physical Volume("cell", 31) = {1};
//+
Physical Volume("strain_gauge", 32) = {2};
//+
Physical Curve("load", 33) = {15};
//+
Physical Surface("fixed", 34) = {4};

// Generar malla al abrir
Mesh 3;