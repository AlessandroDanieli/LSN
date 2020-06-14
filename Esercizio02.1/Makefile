CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS = main.o auxiliary.o random.o stat.o function.o

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
clean:
	rm -rf *.o main.exe *.out
all:
	@make -s
	@./main.exe
