CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS = main.o random.o canonical.o stat.o

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
	
clean:
	rm -rf *.out config.final *.o
	rm -rf OutputData/*
final:
	rm -rf *.o
all:
	@make -s
	@./main.exe
	

