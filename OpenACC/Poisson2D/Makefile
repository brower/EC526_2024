# Copyright (c) 2016, NVIDIA CORPORATION. All rights reserved.

#ACCEL_TYPE=PGI-multicore
ACCEL_TYPE=PGI-tesla
#ACCEL_TYPE=PGI

CC = pgcc
CFLAGS = -DUSE_DOUBLE

ifeq ($(ACCEL_TYPE),PGI)
	CFLAGS += -fast
else ifeq ($(ACCEL_TYPE),PGI-tesla)
	CFLAGS += -Minfo=accel -fast -acc -ta=tesla:cc60
else ifeq ($(ACCEL_TYPE),PGI-multicore)
	CFLAGS += -Minfo=accel -fast -acc -ta=multicore
endif

v0: poisson2d_v0
v1: poisson2d_v1
v2: poisson2d_v2
v3: poisson2d_v3

poisson2d_reference.o: poisson2d_reference.c common.h  Makefile
	$(CC) -c $(CFLAGS) poisson2d_reference.c -o $@

poisson2d_v0: poisson2d_v0.c common.h poisson2d_reference.o Makefile
	$(CC) $(CFLAGS) poisson2d_v0.c poisson2d_reference.o -o $@

poisson2d_v1: poisson2d_v1.c common.h poisson2d_reference.o Makefile
	$(CC) $(CFLAGS) poisson2d_v1.c poisson2d_reference.o -o $@

poisson2d_v2: poisson2d_v2.c common.h poisson2d_reference.o Makefile
	$(CC) $(CFLAGS) poisson2d_v2.c poisson2d_reference.o -o $@

poisson2d_v3: poisson2d_v3.c common.h poisson2d_reference.o Makefile
	$(CC) $(CFLAGS) poisson2d_v3.c poisson2d_reference.o -o $@


.PHONY: clean 

clean:
	rm -f poisson2d_v0 poisson2d_v1 poisson2d_v2 poisson2d_v3 poisson2d*.o
