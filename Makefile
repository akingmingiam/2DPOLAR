

FC := gfortran
FFLAGS := -Wall -Wextra -O2 -std=f2008

SRC_DIR := source
HDR_DIR := header
BUILD_DIR := build

SOURCES := \
    mesh_types.F90 \
    fluid_properties.F90 \
    flow_fields.F90 \
    mesh_ops.F90 \
    field_conversion.F90 \
    boundary_update.F90 \
    initial_ops.F90 \
    flux_types.F90 \
    output_writer.F90 \
    time_step_control.F90 \
    riemann_solver.F90 \
    reconstruction_ops.F90\
    flux_ops.F90 \
    field_update_ops.F90 \
    main.F90

OBJECTS := $(addprefix $(BUILD_DIR)/,$(SOURCES:.F90=.o))

vpath %.F90 $(SRC_DIR) $(HDR_DIR)

.PHONY: all clean

all: main

main: $(BUILD_DIR) $(OBJECTS)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -o $@ $(OBJECTS)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: %.F90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -I$(HDR_DIR) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) main