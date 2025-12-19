# =========================================
# Project configuration
# =========================================
TARGET   := load_cell_beam

DIR_GEO  := geo
DIR_STEP := step
DIR_MSH  := mesh

GEO      := $(DIR_GEO)/$(TARGET).geo
OPT      := $(DIR_GEO)/$(TARGET).opt
STEP     := $(DIR_STEP)/$(TARGET).step
MSH      := $(DIR_MSH)/$(TARGET).msh

# Detect OS (macOS / Linux / Windows)
OS := $(shell uname 2>/dev/null || echo Windows)
PARAVIEW_APP := /Applications/ParaView-6.0.0-RC1.app/Contents/MacOS/paraview


# Gmsh executable (override if needed)
GMSH := gmsh
ifeq ($(OS),Windows)
	GMSH := gmsh.exe
endif

# =========================================
# Targets
# =========================================

.PHONY: gmsh step mesh clean dirs

# Open Gmsh with GEO + OPT (visual config)
gmsh: dirs
	@echo "▶ Opening Gmsh ($(OS)) with geometry..."
	$(GMSH) "$(GEO)" "$(OPT)"

# Open raw STEP (CAD inspection)
step:
	@echo "▶ Opening STEP file ($(OS))..."
	$(GMSH) "$(STEP)"

# Generate mesh without GUI (batch mode)
mesh: dirs
	@echo "▶ Generating mesh..."
	$(GMSH) "$(GEO)" -3 -o "$(MSH)"

# Create required directories
dirs:
	@mkdir -p $(DIR_GEO) $(DIR_STEP) $(DIR_MSH) 2>/dev/null || true

# Remove generated files
clean:
	@echo "▶ Cleaning generated mesh..."
	@rm -f $(DIR_MSH)/*.msh 2>/dev/null || true

precompile:
	julia --project=. -e 'using Pkg; Pkg.precompile()'

# Run Julia simulation
julia: 
	julia --project=. main.jl

paraview:
	open -a ${PARAVIEW_APP} load_cell_results.vtu 

sim: 
	@echo "▶ Running main simulation..."
	$(MAKE) julia 
	
	