CC = g++
CFLAGS = -ansi -Wall -pedantic -O3 $(INCLUDES)
LIBS = 
INCLUDES =
TARGET = mgsolve

SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

.PHONY : clean depend

clean: 
	@/bin/rm -f $(OBJS) 
	@/bin/rm -f $(TARGET)

depend: 
	@makedepend -- $(CFLAGS) -- $(SRC)