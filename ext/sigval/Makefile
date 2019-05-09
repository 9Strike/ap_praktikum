INC = -Iinc

first: sigval

# As C++ program
sigval: main.o sigval.o
	g++ tmp/main.o tmp/sigval.o -o bin/sigval

main.o:
	g++ -g -c src/main.cpp -o tmp/main.o $(INC)

sigval.o:
	g++ -g -c src/sigval.c -o tmp/sigval.o $(INC)

clean:
	rm bin/sigval
	rm tmp/main.o
	rm tmp/sigval.o

# As python extension
py_sigval: py_sigval_build
	cp $(shell find build -name '*.so') tmp/sigval.so
	rm -r build

py_sigval_build:
	python3 setup.py build

py_clean:
	rm tmp/sigval.so
