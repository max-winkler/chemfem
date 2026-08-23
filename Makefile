CPP_INCLUDE = -Iinclude

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
	src/fem/FEFunction.o \
	src/fem/Element.o \
	src/fem/LagrangeElement.o \
	src/fem/FEExpression.o \
	src/fem/LinearForm.o \
	src/fem/BilinearForm.o \
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
	tests/DiffusionTest.o

TESTS_BIN = $(TESTS:.o=)

%.o: %.cpp
	g++ -c $< ${CPP_INCLUDE} ${CPP_FLAGS} -o $@

tests: $(OBJ) $(TESTS)
	$(foreach TEST,$(TESTS_BIN),g++ $(OBJ) $(TEST).o -o $(TEST);)
clean:
	rm -f src/*/*.o tests/*.o
