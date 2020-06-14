CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS = main.o random.o randomwalk.o stat.o vector.o 

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
clean:
	rm *.o main.exe *.out *.xyz
all:
	@make -s
	@./main.exe
