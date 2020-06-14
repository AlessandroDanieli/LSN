CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS = main.o random.o stat.o position.o needle.o auxiliary.o

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
clean:
	rm *.o main.exe *.out
all:
	@make -s
	@./main.exe
