CC	= cc

# put the one you prefer last
# CFLAGS	= -g -O3 -std=c18 -pedantic -Werror -mcpu=power8 -pthread
# CFLAGS	= -g -O3 -std=c18 -pedantic -Werror
# CFLAGS	= -std=c11 -pedantic -Wall -Werror -g -Os -m64 -funsigned-char -static
# CFLAGS	= -std=c11 -pedantic -Wall -g -Os -m64 -funsigned-char -static
CFLAGS		= -g -O3

LDFLAGS		= -g -O3 -pthread -lm

OBJS		= intopt.o main.o error.o

main: clean $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -DPERFORMANCE
	rm -f score
	./a.out -t 4

record:
	sudo perf record -F 10000 -g ./a.out < test/20/0/11637000/i
	sudo perf script | ./FlameGraph/stackcollapse-perf.pl > out.perf-folded
	./FlameGraph/flamegraph.pl out.perf-folded > perf.svg
	google-chrome-stable perf.svg

# Rule to run the test
test: build
	@if [ -z "$(INPUT_DIR)" ]; then \
		echo "ERROR: INPUT_DIR is not set. Usage: make test INPUT_DIR=<path>"; \
		exit 1; \
	fi
	./a.out -t 4 $(INPUT_DIR)/i > sol
	diff sol <(head -n 1 $(INPUT_DIR)/intopt.sol)

# Rule to clean build
clean:
	rm -f *.o a.out

# Rule to build the intopt file
build: clean
	gcc $(CFLAGS) -DINCLUDE_MAIN intopt.c -lm

debug:
	gcc $(CFLAGS) -fsanitize=address -DINCLUDE_MAIN intopt.c -lm
	./a.out

# for use on power.cs.lth.se and other linux machines with oprofile
operf: a.out
	operf -e CYCLES:100000:0:0:1 ./a.out -t 120
