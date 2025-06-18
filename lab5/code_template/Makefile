CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
LDFLAGS = -lm

all: color_transfer

color_transfer: main.cpp stb_image.h stb_image_write.h
	$(CXX) $(CXXFLAGS) -o color_transfer main.cpp $(LDFLAGS)

clean:
	rm -f color_transfer color_transfer_result.png 
