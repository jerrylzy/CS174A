OPTIMIZE = -O2

CC = g++
CFLAGS = $(OPTIMIZE) -std=c++11 -march=native -mtune=native

raytrace: raytrace.cpp
	$(CC) $(CFLAGS) raytrace.cpp -o $@

copy:
	cp raytrace ../../raytrace

clean:
	rm -f *.o *.s raytrace
