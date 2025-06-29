CXX = c++ -std=c++11
CXXFLAGS = -O3 -Wall -Wextra
PINEAPPL_DEPS != pkg-config --cflags --libs pineappl_capi
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf

PROGRAMS = \
	fill-grid-deprecated \
	fill-grid \
	fill-custom-grid-deprecated \
	advanced-convolution-deprecated \
	advanced-convolution \
	advanced-filling-deprecated \
	advanced-filling \
	convolve-grid-deprecated \
	convolve-grid \
	get-subgrids \
	evolve-grid \
	evolve-grid-identity \
	deprecated \
	display-channels-deprecated \
	display-channels \
	display-orders \
	display-orders-deprecated \
	merge-grids \
	modify-grid

all: $(PROGRAMS)

test-examples: $(PROGRAMS)
	set -e && for i in $(PROGRAMS); do diff -u $${i}.output <(./$${i}); done

advanced-convolution-deprecated: advanced-convolution-deprecated.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

advanced-convolution: advanced-convolution.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

advanced-filling-deprecated: advanced-filling-deprecated.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

advanced-filling: advanced-filling.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

convolve-grid-deprecated: convolve-grid-deprecated.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

convolve-grid: convolve-grid.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

deprecated: deprecated.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

get-subgrids: get-subgrids.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

evolve-grid: evolve-grid.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

evolve-grid-identity: evolve-grid-identity.cpp
	$(CXX) $(CXXFLAGS) $< $(LHAPDF_DEPS) $(PINEAPPL_DEPS) -o $@

display-channels-deprecated: display-channels-deprecated.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

display-channels: display-channels.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

display-orders-deprecated: display-orders-deprecated.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

display-orders: display-orders.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

fill-custom-grid-deprecated: fill-grid-deprecated.cpp
	$(CXX) $(CXXFLAGS) $< -DUSE_CUSTOM_GRID_PARAMETERS $(PINEAPPL_DEPS) -o $@

fill-grid-deprecated: fill-grid-deprecated.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

fill-grid: fill-grid.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

merge-grids: merge-grids.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

modify-grid: modify-grid.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) -o $@

PHONY: clean

clean:
	rm -f $(PROGRAMS) *.pineappl *.lz4
