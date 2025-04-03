# Variabili
CC = mpicc #era: gcc
CFLAGS = -Wall -g -lm

# Target predefinito
all: gsa

# Regole per la compilazione
gsa: main.o parallel_gca.o serial_GSA.o utility.o test_functions.o merge_sort.o common.o
	$(CC) $(CFLAGS) -o gsa main.o parallel_gca.o serial_GSA.o utility.o test_functions.o merge_sort.o common.o

# Regole per compilare i file oggetto
main.o: main.c parallel_gca.h serial_GSA.h utility.h test_functions.h merge_sort.h common.h
	$(CC) $(CFLAGS) -c main.c


parallel_gca.o: parallel_gca.c parallel_gca.h
	$(CC) $(CFLAGS) -c parallel_gca.c

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

# Pulizia dei file generati
clean:
	rm -f gsa *.o