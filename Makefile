CPP_INCLUDE = -Iinclude
CPP_FLAGS   = -g

src/linalg/SparseMatrix.o: src/linalg/SparseMatrix.cpp
	g++ -c src/linalg/SparseMatrix.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o src/linalg/SparseMatrix.o

SparseMatrixTest: tests/SparseMatrixTest.cpp src/linalg/SparseMatrix.o
	g++ -c tests/SparseMatrixTest.cpp ${CPP_INCLUDE} ${CPP_FLAGS} -o tests/SparseMatrixTest.o
	g++ tests/SparseMatrixTest.o src/linalg/SparseMatrix.o -o tests/SparseMatrixTest

clean:
	rm src/*.o tests/*.o
