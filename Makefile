CC=g++
CFLAGS=-fopenmp -g -O3 -std=gnu++0x -I./../../boost/ -I./../qft++/include -I ./../eigen/
LDFLAGS=-fopenmp
COVERAGE=
DEBUG=
ifneq ($(COVERAGE),)
  CFLAGS += --coverage -fprofile-arcs -ftest-coverage
  LDFLAGS += -fprofile-arcs
endif
ifneq ($(DEBUG),)
  CFLAGS += -g
endif


C_OBJS=source/DSE/Quark_parameters.cpp
#C_OBJS+=Abs/AbsDiagram.cpp
C_OBJS+=source/DSE/Propagator.cpp
#C_OBJS+=Kernel/Gluon.cpp  $(C_OBJS)

SOURCES=source/main_test.cpp# $(C_OBJS)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=bin/main_test 

SOURCES_T=UnitTests/Integration_test.cpp $(C_OBJS)
OBJECTS_T=$(SOURCES_T:.cpp=.o)
EXECUTABLE_T=bin/Integration_test


all: $(EXECUTABLE) $(EXECUTABLE_T)

	
$(EXECUTABLE_T): $(OBJECTS_T) 
	$(CC) $(LDFLAGS) $(OBJECTS_T) -o $@

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -c -o $@


Integration_Test: $(EXECUTABLE_T)
	./bin/Integration_test --run_test=IntegrationTest
#--report_level=detailed
Unit_Test: $(EXECUTABLE_T)
	./bin/Integration_test --run_test=Gau*

Path_Unit_Test: $(EXECUTABLE_T)
	./bin/Integration_test --run_test=Geom*

clean:
	$(RM) $(EXECUTABLE) $(EXECUTABLE_T)
	$(RM) $(OBJECTS) $(OBJECTS_T)

.PHONY: clean test all
