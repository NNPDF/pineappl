CXX = c++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra
PINEAPPL_DEPS != pkg-config --cflags --libs pineappl_capi
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf

all: convolute-grid fill-custom-grid fill-grid

test-examples: convolute-grid fill-custom-grid fill-grid
	./fill-grid
	./fill-custom-grid
	./convolute-grid

convolute-grid: convolute-grid.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

fill-custom-grid: fill-grid.cpp
	$(CXX) $(CXXFLAGS) $< -DUSE_CUSTOM_GRID_PARAMETERS $(PINEAPPL_DEPS) -o $@

fill-grid: fill-grid.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

PHONY: clean

clean:
	rm -f convolute-grid fill-custom-grid fill-grid
