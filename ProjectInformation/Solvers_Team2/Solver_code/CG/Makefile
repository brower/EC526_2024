.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================
TARGET	=  cg2D

C_SOURCES =  cg2D.cpp
C_OBJS     =  cg2D.o 
MY_INCLUDES =
#triangle.h


CCX = g++
CXXFLAGS = -g -Wall

#============================================================
all: $(TARGET)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET) :   $(C_OBJS)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================

NOTES = MathematicaAffineFlowClean.nb   data/

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES) $(NOTES)

ALL_FILES = $(ALL_SOURCES)  $(NOTES)



clean:
	rm -f $(TARGET) $(C_OBJS) core *.tar *~*~ 

tar: $(ALL_FILES) 
	tar -cvjf $(TARGET).tar $(ALL_FILES)


$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


