# compiler:
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -O2 -I/usr/X11R6/include
CFLAGS_X11 = -L/usr/X11R6/lib -lX11

# the build target executable:
TARGET = cancer_gillespie_simulation

all: $(TARGET)

$(TARGET): src/main.cpp src/phylogeny.o src/universe.o src/cell.o src/cell_type.o src/Phylogeny.hpp src/Universe.hpp src/Cell.hpp src/CellType.hpp
	$(CC) $(CFLAGS) $(CFLAGS_X11) src/phylogeny.o src/universe.o src/cell.o src/cell_type.o src/main.cpp -o $(TARGET)

src/phylogeny.o: src/Phylogeny.cpp src/Phylogeny.hpp
	$(CC) $(CFLAGS) -c src/Phylogeny.cpp -o src/phylogeny.o

src/universe.o: src/Universe.cpp src/Universe.hpp
	$(CC) $(CFLAGS) -c src/Universe.cpp -o src/universe.o

src/cell.o: src/Cell.cpp src/Cell.hpp
	$(CC) $(CFLAGS) -c src/Cell.cpp -o src/cell.o

src/cell_type.o: src/CellType.cpp src/CellType.hpp
	$(CC) $(CFLAGS) -c src/CellType.cpp -o src/cell_type.o

clean:
	$(RM) $(TARGET) src/phylogeny.o src/universe.o src/cell.o src/cell_type.o
