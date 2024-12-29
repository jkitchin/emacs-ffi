all:  gsl-linalg.so

clean:
	rm *.o *.so

##############
# GSL linalg #
##############

gsl-linalg.so: gsl-linalg.o 
	gcc -shared -L/opt/homebrew/lib -lgsl -o gsl-linalg.so gsl-linalg.o

gsl-linalg.o: gsl-linalg.c
	gcc -Wall -I/opt/homebrew/include -fPIC -c gsl-linalg.c

test: gsl-linalg.so tests.el
	emacs -batch -q -l tests.el -f test-linalg
	emacs -batch -q -l tests.el -f test-dgemm

