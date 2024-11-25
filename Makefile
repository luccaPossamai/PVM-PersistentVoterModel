INCLUDES = includes
TARGET = pvmL

$(TARGET): $(TARGET).c
	gcc -I$(INCLUDES) -O3 -o run $(TARGET).c -Wall -lm

clean:
	rm -f $(TARGET)
