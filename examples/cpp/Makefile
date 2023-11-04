CXX = c++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra
PINEAPPL_DEPS != pkg-config --cflags --libs pineappl_capi
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf

PROGRAMS = \
	fill-grid \
	fill-custom-grid \
	convolute-grid \
	display-channels \
	display-orders \

test-examples: $(PROGRAMS)
	for i in $(PROGRAMS); do ./$${i}; done

convolute-grid: convolute-grid.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

display-channels: display-channels.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

display-orders: display-orders.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

fill-custom-grid: fill-grid.cpp
	$(CXX) $(CXXFLAGS) $< -DUSE_CUSTOM_GRID_PARAMETERS $(PINEAPPL_DEPS) -o $@

fill-grid: fill-grid.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

PHONY: clean

clean:
	rm -f $(PROGRAMS)
