CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS = main.o random.o brownian.o stat.o auxiliary.o

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
clean:
	rm *.o *.out main.exe
output:
	rm *.o
all:
	@make -s
	@./main.exe
