CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS = main.o random.o ising.o stat.o

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
	
clean:
	rm -rf *.out config.final *.o
final:
	rm -rf ${LIBS}
all:
	@make -s
	@./main.exe
