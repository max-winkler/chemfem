CPP_INCLUDE = -Iinclude
CPP_FLAGS   = -O3 -Wall

OBJ = 	src/linalg/Vector.o \
	src/linalg/SparseMatrix.o \
	src/linalg/SparseMatrixInserter.o \
	src/linalg/DenseMatrix.o \
	src/linalg/IterativeSolver.o \
	src/mesh/Node.o \
	src/mesh/Cell.o \
	src/mesh/Edge.o \
	src/mesh/Mesh.o \
	src/mesh/UnitSquareMesh.o \
	src/mesh/RefData.o \
	src/mesh/RefDataRegular.o \
	src/quadrature/QuadFormula.o \
	src/fem/FESpace.o \
	src/fem/FEFunction.o \
	src/fem/Element.o \
	src/fem/LagrangeElement.o \
	src/fem/FEExpression.o \
	src/fem/LinearForm.o \
	src/fem/BilinearForm.o \
	src/fem/ErrorNorm.o

TESTS = tests/SparseMatrixTest.o \
	tests/DenseMatrixTest.o \
	tests/MeshTest.o \
	tests/PoissonTest.o \
	tests/MeshRefineTest.o

TESTS_BIN = $(TESTS:.o=)

%.o: %.cpp
	g++ -c $< ${CPP_INCLUDE} ${CPP_FLAGS} -o $@

tests: $(OBJ) $(TESTS)
	$(foreach TEST,$(TESTS_BIN),g++ $(OBJ) $(TEST).o -o $(TEST);)
clean:
	rm -f src/linalg/*.o src/mesh/*.o src/fem/*.o tests/*.o
