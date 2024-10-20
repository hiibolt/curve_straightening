# Compiler and flags
FC = gfortran
FFLAGS = -O2

# Directories
SRC_DIR = src
SUBC_DIR = src/lib
SUBC_DONE_DIR = src/lib_done

# Source files
SRC = $(wildcard $(SRC_DIR)/*.f $(SRC_DIR)/*.for $(SUBC_DIR)/*.f $(SUBC_DONE_DIR)/*.f)

# Output executable
TARGET = geo2

# Build rule
build:
	$(FC) $(FFLAGS) -o $(TARGET) $(SRC) -lblas

# Clean rule
clean:
	rm geo2.prt
	rm geo00.out
	rm geo00.ps
	rm $(TARGET)

# Test rule
run:
	./$(TARGET)