# source files.
SRC = BucketSort_MergeSort_Papi_Paralela.c

OBJ = $(SRC:.cpp=.o)

OUT = sort

# include directories
INCLUDES = -I. 
 
# C compiler flags 
#CCFLAGS = -O0 -Wall
CCFLAGS = -Og -ggdb -g3 -O3 -Wall -std=c99 -fno-omit-frame-pointer -g -I/share/apps/papi/5.4.1/include -L/share/apps/papi/5.4.1/lib
#CCFLAGS = -O3 -Wall -msse4.1 -march=i686

# compiler
CCC = gcc 
#CCC = /opt/intel/Compiler/11.1/073/bin/ia32/icpc 
#CCC = g++-4.5

# libraries
LIBS = -lm -lpapi -fopenmp

.SUFFIXES: .cpp .c 


default: $(OUT)

.cpp.o:
	$(CCC) $(CCFLAGS) $(INCLUDES)  -c $< -o $@

.c.o:
	$(CCC) $(CCFLAGS) $(INCLUDES) -c $< -o $@

$(OUT): $(OBJ)
	$(CCC) -o $(OUT) $(CCFLAGS) $(OBJ) $(LIBS) 

depend:  dep
#
#dep:
#	makedepend -- $(CFLAGS) -- $(INCLUDES) $(SRC)

clean:
	rm -f *.o .a *~ Makefile.bak $(OUT) *.data *.out perfreport *.old
clear:
	rm -f *.o .a *~ Makefile.bak $(OUT)
