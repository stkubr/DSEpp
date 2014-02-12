CC=g++
CFLAGS=-fopenmp -g -O3 -std=gnu++0x -I./../../boost/ -I./../qft++/include -I ./../eigen/
LDFLAGS=-fopenmp
COVERAGE=0
DEBUG=0

ifeq ($(COVERAGE),1)
  CFLAGS +=  -fprofile-arcs -ftest-coverage
  LDFLAGS += -fprofile-arcs
endif

ifeq ($(DEBUG),1)
  CFLAGS += -g
endif

ifeq ($(DEBUG),1)
  CFLAGS += -g
endif

C_OBJS := $(shell find source/* -name '*.cpp')

C_OBJS_T=UnitTests/IntegrationUnitTest.cpp
C_OBJS_T+=UnitTests/GeometryTests/PathsUnitTest.cpp

SOURCES=source/main_test.cpp $(C_OBJS)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=bin/main_test 

SOURCES_T=UnitTests/UnitTestSuites.cpp $(C_OBJS) $(C_OBJS_T)
OBJECTS_T=$(SOURCES_T:.cpp=.o)
EXECUTABLE_T=bin/UnitTestSuites


all: $(EXECUTABLE) $(EXECUTABLE_T)

	
$(EXECUTABLE_T): $(OBJECTS_T) 
	$(CC) $(LDFLAGS) $(OBJECTS_T) -o $@

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $< -c -o $@


Integration_Test_RL: $(EXECUTABLE_T)
	./bin/UnitTestSuites --run_test=IntegrationTest_Rainbow*
	
Integration_Test_PS: $(EXECUTABLE_T)
	./bin/UnitTestSuites --run_test=IntegrationTest_Pseudo*
	
#--report_level=detailed
Unit_Test: $(EXECUTABLE_T)
	./bin/UnitTestSuites --run_test=Gau*

Path_Unit_Test: $(EXECUTABLE_T)
	./bin/UnitTestSuites --run_test=Geom*

clean:
	$(RM) $(EXECUTABLE) $(EXECUTABLE_T)
	$(RM) $(OBJECTS) $(OBJECTS_T)

.PHONY: clean test all
