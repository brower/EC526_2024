.SUFFIXES:
.SUFFIXES: .o .cpp

TARGET	=  findiff

C_SOURCES =  findiff.cpp
C_OBJS     =  findiff.o
MY_INCLUDES =
# stuff.h

CCX = g++ -O2
#CXXFLAGS = -g -std=c++14
# -Wall

#============================================================
all: $(TARGET)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET) :   $(C_OBJS)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@

# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)

NOTES = SeriesExpansion.nb
clean:
	rm -f $(TARGET) $(C_OBJS) core 

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf $(TARGET).tar $(ALL_SOURCES)  $(NOTES)

$(TARGET).ps: $(ALL SOURCES)
	enscript -pcode.ps $(ALL_SOURCES)


