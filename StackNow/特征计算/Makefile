CXX = g++
CXXFLAGS = -O2 -Wall `root-config --cflags`
ROOTLIBS = `root-config --libs` -lMinuit

TARGET = main
OBJS = main.o NKGUtils.o

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(ROOTLIBS)

main.o: main.cpp NKGUtils.h
	$(CXX) $(CXXFLAGS) -c main.cpp

NKGUtils.o: NKGUtils.cpp NKGUtils.h
	$(CXX) $(CXXFLAGS) -c NKGUtils.cpp

clean:
	rm -f *.o $(TARGET)

