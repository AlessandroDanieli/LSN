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
	
#// Jupiter variare lunghezza blocchi e vedere quando "satura", allora numero di blocchi OK
#// Equilibrare sempre prima della simulazione!!!
#// Grandezze scoordinate per montecarlo per un ordine di grandezza ~100 (<~1000)


# Move() already measures u, p (changes)
# Add SampleG()
# Tail corrections?
# Acceptance ratio 
