# Paths that differ from machine to machine. Make.inc is not part of the repository:
# copy Make.inc.example to Make.inc and adjust it where the defaults below do not fit.
# The leading dash keeps make quiet when the file is absent.
-include Make.inc

# Defaults for a distribution install of SuiteSparse. Assignments made in Make.inc
# take precedence, also when they set a variable to the empty string.
UMFPACK_INCLUDE ?= -I/usr/include/suitesparse
UMFPACK_LIB     ?= -lumfpack

CPP_INCLUDE = -Iinclude ${UMFPACK_INCLUDE}

CPP_DEBUG_FLAGS   = -g -O0 -Wall -Wextra -fno-omit-frame-pointer

# -march=native ties the binaries to the CPU they were built on, and -ffast-math
# reorders floating point arithmetic and drops the handling of NaN and Inf
CPP_RELEASE_FLAGS = -O3 -march=native -mtune=native -ffast-math -funroll-loops \
                    -DNDEBUG -Wall

CPP_FLAGS = ${CPP_RELEASE_FLAGS}

OBJ = 	src/linalg/Vector.o \
	src/linalg/SparseMatrix.o \
	src/linalg/SparseMatrixInserter.o \
	src/linalg/DenseMatrix.o \
	src/linalg/IterativeSolver.o \
	src/linalg/DirectSolver.o \
	src/mesh/Node.o \
	src/mesh/Cell.o \
	src/mesh/Edge.o \
	src/mesh//CellInfo.o \
	src/mesh/Mesh.o \
	src/mesh/UnitSquareMesh.o \
	src/mesh/LShapeMesh.o \
	src/mesh/RefData.o \
	src/mesh/RefDataRegular.o \
	src/mesh/RefDataBisection0.o \
	src/mesh/RefDataBisection1.o \
	src/mesh/RefDataBisection2.o \
	src/quadrature/QuadFormula.o \
	src/fem/FESpace.o \
	src/fem/DirichletValues.o \
	src/fem/DofManager.o \
	src/fem/DofTransform.o \
	src/fem/FEFunction.o \
	src/fem/Element.o \
	src/fem/LagrangeElement.o \
	src/fem/HermiteElement.o \
	src/fem/ProductElement.o \
	src/fem/CrouzeixRaviartElement.o \
	src/fem/DGElement.o \
	src/fem/RaviartThomasElement.o \
	src/fem/FEExpression.o \
	src/fem/LinearForm.o \
	src/fem/BilinearForm.o \
	src/fem/PointValues.o \
	src/fem/BlockSystem.o \
	src/fem/VtkOutput.o \
	src/fem/ErrorNorm.o \
	src/fem/GenericEstimator.o	

TESTS = tests/SparseMatrixTest.o \
	tests/DenseMatrixTest.o \
	tests/MeshTest.o \
	tests/PoissonTest.o \
	tests/MeshRefineTest.o \
	tests/LocalRefinement.o \
	tests/AdaptivityTest.o \
	tests/ConvectionTest.o \
	tests/EstimatorTest.o \
	tests/DiffusionTest.o \
	tests/SolverTest.o \
	tests/NeumannTest.o \
	tests/CrouzeixRaviartTest.o \
	tests/StokesTest.o \
	tests/InstationaryStokesTest.o \
	tests/RobinTest.o \
	tests/VectorPoissonTest.o \
	tests/DGTest.o \
	tests/RaviartThomasTest.o \
	tests/MixedPoissonTest.o \
	tests/HermiteTest.o

TESTS_BIN = $(TESTS:.o=)

%.o: %.cpp
	g++ -c $< ${CPP_INCLUDE} ${CPP_FLAGS} -o $@

tests: $(OBJ) $(TESTS)
	$(foreach TEST,$(TESTS_BIN),g++ $(OBJ) $(TEST).o -o $(TEST) ${UMFPACK_LIB} &&) true

# Object files do not record the flags they were built with, so switching the
# build type has to start from scratch
debug:
	$(MAKE) clean
	$(MAKE) CPP_FLAGS="${CPP_DEBUG_FLAGS}" tests

release:
	$(MAKE) clean
	$(MAKE) CPP_FLAGS="${CPP_RELEASE_FLAGS}" tests

clean:
	rm -f src/*/*.o tests/*.o

.PHONY: debug release clean
