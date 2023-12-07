CC	= cc
CFLAGS	= -std=c11 -pedantic -Wall -Werror -g -O3 -funsigned-char
OUT	= a.out
FILE	= poly

execute: compile
	./a.out > output
	diff output correct

test: compile
	./a.out

own: clean $(FILE).o
	$(CC) $(CFLAGS) -DINCLUDE_MAIN -o $(OUT) $(FILE).c
	./a.out

compile: clean $(FILE).o malloc.o main.o
	$(CC) $(CFLAGS) -o $(OUT) $(FILE).o malloc.o main.o

$(FILE).o: $(FILE).c
	$(CC) -include wecc.h -c $(FILE).c $(CFLAGS)

malloc.o: malloc.c
	$(CC) -c malloc.c $(CFLAGS) 

clean:
	rm -f *.s *.o $(OUT) output
