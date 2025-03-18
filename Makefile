# Variabili
CC = mpicc #era: gcc
CFLAGS = -Wall -g -lm

# Target predefinito
all: gca

# Regole per la compilazione
gca: main.o parallel_gca.o serial_gca.o utility.o test_functions.o merge_sort.o
	$(CC) $(CFLAGS) -o gca main.o parallel_gca.o serial_gca.o utility.o test_functions.o merge_sort.o

# Regole per compilare i file oggetto
main.o: main.c parallel_gca.h serial_gca.h
	$(CC) $(CFLAGS) -c main.c

parallel_gca.o: parallel_gca.c parallel_gca.h
	$(CC) $(CFLAGS) -c parallel_gca.c

serial_gca.o: serial_gca.c serial_gca.h
	$(CC) $(CFLAGS) -c serial_gca.c

utility.o: utility.c utility.h
	$(CC) $(CFLAGS) -c utility.c

test_functions.o: test_functions.c test_functions.h
	$(CC) $(CFLAGS) -c test_functions.c

merge_sort.o: merge_sort.c merge_sort.h
	$(CC) $(CFLAGS) -c merge_sort.c

# Pulizia dei file generati
clean:
	rm -f gca *.o