INCLUDES = headers
LIB = lib
TARGET = pvmL

comp: $(TARGET).c
	gcc -I$(INCLUDES) -L$(LIB) -O3 -o run $(TARGET).c -Wall -llplibs -lm
comprun: $(TARGET).c
	gcc -I$(INCLUDES) -L$(LIB) -O3 -o run $(TARGET).c -Wall -llplibs -lm
	./run
debug: $(TARGET).c
	gcc -FORCE_SEED=1 -g -I$(INCLUDES) -L$(LIB) -O3 -o run $(TARGET).c -Wall -llplibs -lm
	./run
plot: $(TARGET).c
	gcc -DFORCE_SEED=123456789 -DPLOT=1 -I$(INCLUDES) -L$(LIB) -O3 -o run $(TARGET).c -Wall -llplibs -lm
	./run | gnuplot
clean:
	rm -f run
