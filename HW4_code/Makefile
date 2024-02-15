.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================



TARGET2 = gaussQuad
C_OBJS2 = gaussQuadWeight.o

TARGET3 = multiDim
C_OBJS3 = multiDimIntegrals.o


ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)

DATA_FILES =
#x8err.pdf  cosPIxerr.pdf x2p1inverr.pdf 

CCX = g++
CXXFLAGS = -g -Wall
#-std=c99

#============================================================
all: $(TARGET1) $(TARGET2) $(TARGET3)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET1) :   $(C_OBJS1)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@
$(TARGET2) :   $(C_OBJS2)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@ 
$(TARGET3) :   $(C_OBJS3)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@ 


# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================


clean:
	rm -f $(TARGET1) $(TARGET2) $(TARGET3)  $(C_OBJS1) $(C_OBJS2)  $(C_OBJS3) *~

tar: $(ALL_SOURCES) $(DATA_FILES)
	tar cvf HW3_code.tar $(ALL_SOURCES)  $(DATA_FILES)



