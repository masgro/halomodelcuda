VPATH := src/

EXEC   = sigma

OBJS   = main.o power.o allvars.o read_param.o deltas.o polifitgsl.o

INCL   = allvars.h proto.h Makefile

CC       =   gcc        # sets the C-compiler (default)
OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)
MPICHLIB = 

LIBS   =   -lm  -lgsl -lgslcblas

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



