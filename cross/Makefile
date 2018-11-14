# Binary file
BIN = hm.cono

OPTIONS := -DCG
OPTIONS += -DMERCHAN
#OPTIONS += -DNEW_ALIGN
OPTIONS += -DTAO_II

ifeq ($(term),2h)
  OPTIONS += -DDOSH
	BIN1 = hm.cono.2h.x
	BIN2 = reescribe.2h.x
else
	BIN1 = hm.cono.x
	BIN2 = reescribe.x
endif

#OPTIONS += -DCENTROS_MASA_MIN=12.23
#OPTIONS += -DCENTROS_MASA_MAX=13.72

#OPTIONS += -DCENTROS_MASA_MIN=13.72
#OPTIONS += -DCENTROS_MASA_MAX=14.14

#OPTIONS += -DCENTROS_MASA_MIN=14.14
#OPTIONS += -DCENTROS_MASA_MAX=15.34

OPTIONS += -DANGULO=45
OPTIONS += -DMAG=19

# Flags
CFLAGS = -O3 --use_fast_math #--ptxas-options=-v
LDFLAGS = -lcudart -lgomp -lgsl -lgslcblas -lm

MACHINE=$(shell hostname)

# Compilers
CC = nvcc
LINKER = nvcc

ifeq ($(MACHINE),mendieta)
OPTIONS += -DDEVICE=0
LDFLAGS  = -lcuda -lcudart -lgomp -lgsl -lgslcblas -lm
CFLAGS  += -arch=sm_20
endif

ifeq ($(MACHINE),clemente)
OPTIONS += -DDEVICE=0
LDFLAGS  = -lcuda -lcudart -lgomp -lgsl -lgslcblas -lm
CFLAGS  += -arch=sm_20
endif

# Files
MAKEFILE = Makefile
C_SOURCES = hm.cono.cu
HEADERS =
C_OBJS = $(patsubst %.cu, %.o,$(C_SOURCES))

todo: $(BIN) reescribe.x

# Rules
$(BIN): $(C_OBJS) $(HEADERS) $(MAKEFILE)
	$(LINKER) $(C_OBJS) $(LDFLAGS) $(INCLUDES) $(LIBS) -o $(BIN1)

$(C_OBJS): $(C_SOURCES) $(HEADERS) $(MAKEFILE)
	$(CC) $(OPTIONS) -c $(C_SOURCES) $(CFLAGS) $(INCLUDES)

ngmedio.x: ngmedio.cu $(HEADERS) $(MAKEFILE)
	$(LINKER) $(OPTIONS) $< $(LDFLAGS) $(INCLUDES) $(LIBS) -o $@

reescribe.o: reescribe.c $(MAKEFILE)
	gcc $(OPTIONS) reescribe.c -c

reescribe.x: reescribe.o $(MAKEFILE)
	gcc $(OPTIONS) reescribe.o -lgsl -lgslcblas -lm -o $(BIN2)

clean:
	rm -f *.o
