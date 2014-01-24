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
C_OBJS+=source/DedicMem/DedicMem.cpp
#C_OBJS+=Abs/AbsDiagram.cpp
C_OBJS+=source/DSE/Propagator.cpp
C_OBJS+=source/DSE/Quark.cpp
C_OBJS+=source/NumLibs/Geometry/Path.cpp
C_OBJS+=source/NumLibs/Geometry/Line.cpp
C_OBJS+=source/NumLibs/Geometry/Parabola.cpp
C_OBJS+=source/NumLibs/Geometry/ParabolaContour.cpp
#C_OBJS+=Kernel/Gluon.cpp  $(C_OBJS)

C_OBJS_T=UnitTests/IntegrationUnitTest.cpp
C_OBJS_T+=UnitTests/GeometryTests/PathsUnitTest.cpp

SOURCES=source/main_test.cpp $(C_OBJS)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=bin/main_test 

SOURCES_T=$(C_OBJS) UnitTests/UnitTestSuites.cpp $(C_OBJS_T)
OBJECTS_T=$(SOURCES_T:.cpp=.o)
EXECUTABLE_T=bin/UnitTestSuites


all: $(EXECUTABLE) $(EXECUTABLE_T)

	
$(EXECUTABLE_T): $(OBJECTS_T) 
	$(CC) $(LDFLAGS) $(OBJECTS_T) -o $@

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -c -o $@


Integration_Test: $(EXECUTABLE_T)
	./bin/UnitTestSuites --run_test=IntegrationTest
#--report_level=detailed
Unit_Test: $(EXECUTABLE_T)
	./bin/UnitTestSuites --run_test=Gau*

Path_Unit_Test: $(EXECUTABLE_T)
	./bin/UnitTestSuites --run_test=Geom*

clean:
	$(RM) $(EXECUTABLE) $(EXECUTABLE_T)
	$(RM) $(OBJECTS) $(OBJECTS_T)

.PHONY: clean test all
