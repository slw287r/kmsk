ARCH := $(shell arch)
ifeq ($(ARCH),x86_64)
CFLAGS=-Wall -O3 -static -lz -lpthread
else
CFLAGS=-Wall -O3 -lz -lpthread
endif

kmsk: kmsk.c cgranges.c thpool.c
	cc -o $@ $^ $(CFLAGS)

clean:
	rm -f kmsk
