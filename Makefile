CPP_INCLUDE = -Iinclude
CPP_FLAGS   = -g

OBJ = 	src/linalg/Vector.o \
	src/linalg/SparseMatrix.o \
	src/linalg/SparseMatrixInserter.o \
	src/mesh/Node.o \
	src/mesh/Cell.o \
	src/mesh/Mesh.o \
	src/mesh/UnitSquareMesh.o \

%.o: %.cpp
	g++ -c $< ${CPP_INCLUDE} ${CPP_FLAGS} -o $@

src/linalg/SparseMatrixInserter.o: src/linalg/SparseMatrixInserter.cpp
	g++ -c src/linalg/SparseMatrixInserter.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o src/linalg/SparseMatrixInserter.o

tests: tests/SparseMatrixTest.cpp tests/MeshTest.cpp $(OBJ)
	g++ -c tests/MeshTest.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o tests/MeshTest.o
	g++ -c tests/SparseMatrixTest.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o tests/SparseMatrixTest.o
	g++ tests/SparseMatrixTest.o $(OBJ) -o tests/SparseMatrixTest
	g++ tests/MeshTest.o $(OBJ) -o tests/MeshTest

clean:
	rm src/linalg/*.o src/mesh/*.o tests/*.o
