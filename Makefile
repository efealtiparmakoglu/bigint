CC      ?= cc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -Wpedantic

all: bigint

bigint: src/main.c src/bigint.c src/bigint.h
	$(CC) $(CFLAGS) -o $@ src/main.c src/bigint.c

bigint_tests: tests/bigint_tests.c src/bigint.c src/bigint.h
	$(CC) $(CFLAGS) -Isrc -o bigint_tests tests/bigint_tests.c src/bigint.c

test: bigint_tests bigint
	./bigint_tests
	./bigint fib 300 | grep -q '^222232244629420445529739893461909967206666939096499764990979600' && echo "fib OK"
	./bigint fact 40 | grep -q '^815915283247897734345611269596115894272000000000' && echo "fact OK"

clean:
	rm -f bigint bigint_tests

.PHONY: all test clean
