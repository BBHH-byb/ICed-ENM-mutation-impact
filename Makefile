FC = gfortran
OPENBLAS_PREFIX ?= /opt/homebrew/opt/openblas

FFLAGS ?= -O3 -ffree-line-length-none
OMPFLAGS ?= -fopenmp
BLAS_LIBS ?= -L$(OPENBLAS_PREFIX)/lib -lopenblas -Wl,-rpath,$(OPENBLAS_PREFIX)/lib

SRC_DIR := src
BIN_DIR := bin
TARGET := $(BIN_DIR)/ICed_ENM_NMA

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC_DIR)/ICed_ENM_NMA.f95
	mkdir -p $(BIN_DIR)
	$(FC) $(OMPFLAGS) $(FFLAGS) $< -o $@ $(BLAS_LIBS)

clean:
	rm -rf $(BIN_DIR)
