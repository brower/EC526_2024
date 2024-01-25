.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================
TARGET	= diff

C_SOURCES =  findiff_soln.cpp
C_OBJS     =  findiff_soln.o 
MY_INCLUDES = 


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

NOTES =


ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES) $(NOTES)

ALL_FILES = $(ALL_SOURCES)  $(NOTES)

#Jacobian.nb


clean:
	rm -f $(TARGET) $(C_OBJS) core *.tar *~*~ 

tar: $(ALL_FILES) 
	tar -cvjf $(TARGET).tar $(ALL_FILES)


$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


