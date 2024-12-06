INCLUDES = includes
TARGET = pvmL

comp: $(TARGET).c
	gcc -I$(INCLUDES) -O3 -o run $(TARGET).c -Wall -lm
	./run
debug: $(TARGET).c
	gcc -FORCE_SEED=1 -I$(INCLUDES) -O3 -o run $(TARGET).c -Wall -lm
	./run
clean:
	rm -f $(TARGET)
