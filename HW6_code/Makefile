.SUFFIXES:
.SUFFIXES: .o .cpp

#============================================================
TARGET	=  mg

C_SOURCES = mg.cpp
C_OBJS     =  mg.o
MY_INCLUDES =


CCX = g++
CXXFLAGS = -g -O2  $(INC)

#============================================================
all: $(TARGET)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET) :   $(C_OBJS)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites 
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)

NOTES = FFT_2021.pdf

clean:
	rm -f $(TARGET) $(C_OBJS) *.*~ core 

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf $(TARGET).tar $(ALL_SOURCES)  $(NOTES)


