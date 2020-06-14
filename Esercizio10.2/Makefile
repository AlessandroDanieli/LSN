CC= mpicxx
CFLAGS = -Wall -O3 --std=c++11
LIBS = main.o random.o genetic.o vector.o

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
	
clean:
	rm -rf *.out config.final *.o
final:
	rm -rf *.o
all:
	@make -s
	@mpirun -np 4 ./main.exe
