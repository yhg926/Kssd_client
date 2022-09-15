CC=gcc
CFLAGS= -std=gnu11 -Wall -Wno-unused-result -O3 -ggdb -lz -fopenmp
all:
	$(CC) $(CFLAGS)  *.c -o ./kssd_client -lm
alert:
	$(CC) $(CFLAGS) -DCOMPONENT_SZ=8 *.c -o ./kssd_CSZ8 -lm
