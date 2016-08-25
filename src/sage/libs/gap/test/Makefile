

GAPDIR="$(SAGE_LOCAL)/gap/latest"

# VALGRIND=valgrind --leak-check=full
# VALGRIND=valgrind --db-attach=yes
#VALGRIND=valgrind --suppressions=libgap.supp --gen-suppressions=yes
#VALGRIND=valgrind --suppressions=libgap.supp --db-attach=yes
VALGRIND=

all: main
	$(VALGRIND) ./main

main: main.o Makefile
	gcc -L$(SAGE_LOCAL)/lib -o main main.o -lgap -lcsage -lntl -lstdc++ -lpari -lpython2.7 -lm -lgmp

main.o: main.c
	echo $(GAPDIR)
	gcc -std=gnu99 -DGAPDIR="\"$(GAPDIR)\"" -I$(SAGE_LOCAL)/include -I$(SAGE_LOCAL)/include/python2.6 -c -g $^

clean:
	rm main.o main *~

.PHONY: all clean
