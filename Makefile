CC=gcc
CFLAGS=-g -fsanitize=address -fopenmp -std=c11 -Wall -Werror -O0
TARGET=pagerank
.PHONY: clean
all: $(TARGET)

pagerank: ./src/pagerank.c ./src/lib/pagerank.h
	$(CC) $(CFLAGS) $^ -o $@ -lpthread -lm

benchmark:
	make clean
	make pagerank
	time bash tester.sh
	make clean

clean:
	rm -f *.o
	rm -f pagerank
