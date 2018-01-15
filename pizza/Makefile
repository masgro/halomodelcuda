# Binary file
BIN = hm.cono

OPTIONS := $(kind)
OPTIONS += -DMERCHAN
OPTIONS += -DTAO_II
#OPTIONS += -DSLOAN
OPTIONS += -DNUM=$(num)

ifeq ($(term),2h)
  OPTIONS += -DDOSH
endif

ifndef bcmedio
OPTIONS += -DBCMEDIO=0.8
else
OPTIONS += -DBCMEDIO=$(bcmedio)
endif

ifndef abmedio
OPTIONS += -DABMEDIO=0.7
else
OPTIONS += -DABMEDIO=$(abmedio)
endif

ifndef align_b
OPTIONS += -DALIGN_B=0.3
else
OPTIONS += -DALIGN_B=$(align_b)
endif

ifndef align_c
OPTIONS += -DALIGN_C=0.5
else
OPTIONS += -DALIGN_C=$(align_c)
endif

ifndef m1
OPTIONS += -DCENTROS_MASA_MIN=12.0
else
OPTIONS += -DCENTROS_MASA_MIN=$(m1)
endif

ifndef m2
OPTIONS += -DCENTROS_MASA_MAX=12.5
else
OPTIONS += -DCENTROS_MASA_MAX=$(m2)
endif

ifndef angulo
OPTIONS += -DANGULO=45
else
OPTIONS += -DANGULO=$(angulo)
endif

ifndef mag
OPTIONS += -DMAG=19
else
OPTIONS += -DMAG=$(mag)
endif

# Flags
CFLAGS = -O3 --use_fast_math #--ptxas-options=-v
LDFLAGS = -lcudart -lgomp -lgsl -lgslcblas -lm

MACHINE=$(shell hostname)
CFLAGS += -arch=sm_20

# Compilers
CC = nvcc
LINKER = nvcc

# Files
MAKEFILE = Makefile
C_SOURCES = $(BIN).cu
HEADERS =
C_OBJS = $(patsubst %.cu, %.o,$(C_SOURCES))

todo: $(BIN) reescribe.x

# Rules
$(BIN): $(C_OBJS) $(HEADERS) $(MAKEFILE)
	$(LINKER) $(C_OBJS) $(LDFLAGS) $(INCLUDES) $(LIBS) -o $(BIN)

$(C_OBJS): $(C_SOURCES) $(HEADERS) $(MAKEFILE)
	$(CC) $(OPTIONS) -c $(C_SOURCES) $(CFLAGS) $(INCLUDES)

reescribe.x: reescribe.c $(MAKEFILE)
	gcc $(OPTIONS) -o $@ reescribe.c -lgsl -lgslcblas -lm

HandleError.o: HandleError.cu
	nvcc -c $< $(CFLAGS)

funciones.o: funciones.cu
	nvcc -c $< $(CFLAGS)

clean:
	rm -f $(BIN) reescribe.x *.o
