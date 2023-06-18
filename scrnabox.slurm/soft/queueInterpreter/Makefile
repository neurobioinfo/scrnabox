CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		QueueInterpreter.o

LIBS =   -lboost_regex -lboost_program_options  -I /usr/include/boost 

TARGET =	QueueInterpreter

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
