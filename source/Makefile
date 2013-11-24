CC=g++
CFLAGS=-fopenmp -O3 -std=gnu++0x -I./../../boost/ -I./../qft++/include -I ./../../eigen/
LDFLAGS=-fopenmp
COVERAGE=
ifneq ($(COVERAGE),)
  CFLAGS += --coverage -fprofile-arcs -ftest-coverage
  LDFLAGS += -fprofile-arcs
endif

#C_OBJS=DSE/Quark_parameters.cpp
#C_OBJS+=Abs/AbsDiagram.cpp
#C_OBJS+=DSE/Propagator.cpp
#C_OBJS+=Kernel/Gluon.cpp  $(C_OBJS)

SOURCES=main_test.cpp $(C_OBJS)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main_test 

SOURCES_T=Integration_test.cpp $(C_OBJS)
OBJECTS_T=$(SOURCES_T:.cpp=.o)
EXECUTABLE_T=Integration_test


all: $(EXECUTABLE) $(EXECUTABLE_T)

	
$(EXECUTABLE_T): $(OBJECTS_T) 
	$(CC) $(LDFLAGS) $(OBJECTS_T) -o $@

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -c -o $@


Integration_Test: $(EXECUTABLE_T)
	./Integration_test --run_test=IntegrationTest
#--report_level=detailed
Unit_Test: $(EXECUTABLE_T)
	./Integration_test --run_test=Gau*

Path_Unit_Test: $(EXECUTABLE_T)
	./Integration_test --run_test=Geom*

clean:
	$(RM) $(EXECUTABLE) $(EXECUTABLE_T)
	$(RM) $(OBJECTS) $(OBJECTS_T)

.PHONY: clean test all
