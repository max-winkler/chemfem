CPP_INCLUDE = -Iinclude
CPP_FLAGS   = -g

OBJ = 	src/linalg/Vector.o \
	src/linalg/SparseMatrix.o \
	src/linalg/SparseMatrixInserter.o \
	src/linalg/DenseMatrix.o \
	src/mesh/Node.o \
	src/mesh/Cell.o \
	src/mesh/Mesh.o \
	src/mesh/UnitSquareMesh.o \
	src/quadrature/QuadFormula.o \
	src/fem/FESpace.o \
	src/fem/Element.o \
	src/fem/LagrangeElement.o \
	src/fem/FEExpression.o \
	src/fem/BilinearForm.o 

TESTS = tests/SparseMatrixTest.o \
	tests/DenseMatrixTest.o \
	tests/MeshTest.o \
	tests/PoissonTest.o

TESTS_BIN = $(TESTS:.o=)

%.o: %.cpp
	g++ -c $< ${CPP_INCLUDE} ${CPP_FLAGS} -o $@


tests: $(TESTS) $(TESTS_BIN)

tests: tests/SparseMatrixTest.cpp tests/MeshTest.cpp \
	tests/DenseMatrixTest.cpp tests/PoissonTest.cpp $(OBJ)
	g++ -c tests/SparseMatrixTest.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o tests/SparseMatrixTest.o
	g++ -c tests/DenseMatrixTest.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o tests/DenseMatrixTest.o
	g++ -c tests/MeshTest.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o tests/MeshTest.o
	g++ -c tests/PoissonTest.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o tests/PoissonTest.o
	g++ tests/SparseMatrixTest.o $(OBJ) -o tests/SparseMatrixTest
	g++ tests/DenseMatrixTest.o $(OBJ) -o tests/DenseMatrixTest
	g++ tests/MeshTest.o $(OBJ) -o tests/MeshTest
	g++ tests/PoissonTest.o $(OBJ) -o tests/PoissonTest

clean:
	rm -f src/linalg/*.o src/mesh/*.o tests/*.o
