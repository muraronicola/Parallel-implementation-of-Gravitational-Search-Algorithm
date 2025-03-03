# Variabili
CC = gcc
CFLAGS = -Wall -g

# Target predefinito
all: gca

# Regole per la compilazione
gca: main.o functions.o
	$(CC) $(CFLAGS) -o gca main.o functions.o

# Regole per compilare i file oggetto
main.o: main.c functions.h
	$(CC) $(CFLAGS) -c main.c

functions.o: functions.c functions.h
	$(CC) $(CFLAGS) -c functions.c

# Pulizia dei file generati
clean:
	rm -f gca *.o