CC = mpicc 
CFLAGS = -Wall -g -lm

# Target 
all: gsa

# Rule for building the executable
gsa: main.o parallel_GSA.o serial_GSA.o utility.o test_functions.o merge_sort.o common.o min_heap.o
	$(CC) $(CFLAGS) -o gsa main.o parallel_GSA.o serial_GSA.o utility.o test_functions.o merge_sort.o common.o min_heap.o

# Rule for compiling object files
main.o: main.c parallel_GSA.h serial_GSA.h utility.h test_functions.h merge_sort.h common.h min_heap.h
	$(CC) $(CFLAGS) -c main.c

parallel_GSA.o: parallel_GSA.c parallel_GSA.h
	$(CC) $(CFLAGS) -c parallel_GSA.c

serial_GSA.o: serial_GSA.c serial_GSA.h
	$(CC) $(CFLAGS) -c serial_GSA.c

utility.o: utility.c utility.h
	$(CC) $(CFLAGS) -c utility.c

test_functions.o: test_functions.c test_functions.h
	$(CC) $(CFLAGS) -c test_functions.c

merge_sort.o: merge_sort.c merge_sort.h
	$(CC) $(CFLAGS) -c merge_sort.c

common.o: common.c common.h
	$(CC) $(CFLAGS) -c common.c


min_heap.c: min_heap.h
	$(CC) $(CFLAGS) -c min_heap.c

# Clean target
clean:
	rm -f gsa *.o