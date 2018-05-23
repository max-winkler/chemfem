CPP_INCLUDE = -Iinclude
CPP_FLAGS   = -g

OBJ = 	src/linalg/SparseMatrix.o	\
	src/linalg/SparseMatrixInserter.o

%.o: %.cpp
	g++ -c $< ${CPP_INCLUDE} ${CPP_FLAGS} -o $@

src/linalg/SparseMatrixInserter.o: src/linalg/SparseMatrixInserter.cpp
	g++ -c src/linalg/SparseMatrixInserter.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o src/linalg/SparseMatrixInserter.o

SparseMatrixTest: tests/SparseMatrixTest.cpp $(OBJ)
	g++ -c tests/SparseMatrixTest.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o tests/SparseMatrixTest.o
	g++ tests/SparseMatrixTest.o $(OBJ) -o tests/SparseMatrixTest

clean:
	rm src/*.o tests/*.o
