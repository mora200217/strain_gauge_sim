SetFactory("OpenCASCADE");

// Importar geometría
Merge "../step/load_cell_beam.step";

// -----------------------------
// Unir volúmenes (CRÍTICO)
// -----------------------------

// Supongamos:
// Volumen 1 → beam
// Volumen 2 → strain gauge

BooleanFragments{ Volume{1}; Delete; }
               { Volume{2}; Delete; };

// -----------------------------
// Physical groups (DESPUÉS del boolean)
// -----------------------------

Physical Volume("cell")         = {1};
Physical Volume("strain_gauge") = {2};

//+
Physical Surface("fixed", 49) = {20};
//+
Physical Surface("load", 50) = {18};
//+
Physical Surface("sensor", 51) = {13};
//+
Physical Surface("sensor_surface", 52) = {14};


// -----------------------------
// Mesh
// -----------------------------

Mesh.ElementOrder = 3;
Mesh.HighOrderOptimize = 2;
Mesh 3;

// Guardar
Save "../mesh/load_cell_beam.msh";


