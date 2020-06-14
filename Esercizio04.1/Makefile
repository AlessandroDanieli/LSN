CC = g++
CFLAGS = -Wall -O3 --std=c++11

LIBS = main.o moldyn.o stat.o
# $@  ->  target name
# $^  ->  

main.exe: ${LIBS}
	$(CC) -o $@ $^ ${CFLAGS}
%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}
clean:
	rm -rf *.o old.0 old.final main.exe
	rm -rf OutputData/output*.dat
	rm -rf frames/*.xyz
	
output:
	rm *.o
all:
	@make -s
	@./main.exe
