a.out : main3.o M3.o ex3.o get_time.o
	g++ main3.o M3.o ex3.o get_time.o -pthread -o a.out -O3 

main3.o: main3.cpp
	g++ -c main3.cpp -o main3.o -O3

M3.o: M3.cpp
	g++ -c M3.cpp -o M3.o -O3
ex3.o: ex3.cpp
	g++ -c ex3.cpp -o ex3.o -O3 -Wno-error=aggressive-loop-optimizations
get_time.o:
	g++ -c get_time.cpp -o get_time.o -O3
clean:
	rm *.o