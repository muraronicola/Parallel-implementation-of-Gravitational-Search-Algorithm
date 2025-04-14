# Compilatore e flag
CC = mpicc
CFLAGS = -Wall -g -lm \
         -Iutil -Imain -Imain/parallel -Imain/serial -Iutil/heap -Iutil/sort

all: gsa

# Rule for building the executable
gsa: main.o main/parallel/parallel_GSA.o main/serial/serial_GSA.o main/common.o util/utility.o util/test_functions.o util/sort/merge_sort.o util/heap/min_heap.o
	$(CC) $(CFLAGS) -o $@ main.o main/parallel/parallel_GSA.o main/serial/serial_GSA.o main/common.o util/utility.o util/test_functions.o util/sort/merge_sort.o util/heap/min_heap.o


# Rule for compiling object files
main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

main/parallel/parallel_GSA.o: main/parallel/parallel_GSA.c main/parallel/parallel_GSA.h
	$(CC) $(CFLAGS) -c main/parallel/parallel_GSA.c -o main/parallel/parallel_GSA.o

main/serial/serial_GSA.o: main/serial/serial_GSA.c main/serial/serial_GSA.h
	$(CC) $(CFLAGS) -c main/serial/serial_GSA.c -o main/serial/serial_GSA.o

main/common.o: main/common.c main/common.h
	$(CC) $(CFLAGS) -c main/common.c -o main/common.o

util/utility.o: util/utility.c util/utility.h
	$(CC) $(CFLAGS) -c util/utility.c -o util/utility.o

util/test_functions.o: util/test_functions.c util/test_functions.h
	$(CC) $(CFLAGS) -c util/test_functions.c -o util/test_functions.o

util/sort/merge_sort.o: util/sort/merge_sort.c util/sort/merge_sort.h
	$(CC) $(CFLAGS) -c util/sort/merge_sort.c -o util/sort/merge_sort.o

util/heap/min_heap.o: util/heap/min_heap.c util/heap/min_heap.h
	$(CC) $(CFLAGS) -c util/heap/min_heap.c -o util/heap/min_heap.o


#Clean target
clean:
	rm -f gsa *.o
