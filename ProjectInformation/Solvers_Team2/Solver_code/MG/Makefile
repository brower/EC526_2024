.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================
TARGET1 =  mg1D
C_OBJS1 = mg1D.o

TARGET2 =  mg2D
C_OBJS2 = mg2D.o


C_SOURCES = mg1D.cpp mg2D.cpp 
C_OBJS     =  mainFFT.o fft.o

MY_INCLUDES =


CCX = g++
CXXFLAGS = -g -O2  $(INC)
#-std=c++11x

#============================================================
all: $(TARGET1) $(TARGET2) 

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<

$(TARGET1) :   $(C_OBJS1)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@
$(TARGET2) :   $(C_OBJS2)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@


# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites 
#============================================================

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)

NOTES = my-plot.pdf

clean:
	rm -f $(TARGET1) $(TARGET2) $(C_OBJS1)  $(C_OBJS2)  *.~*~ core 

tar: $(ALL_SOURCES) $(NOTES)
	tar cvf MG.tar $(ALL_SOURCES)  $(NOTES)


