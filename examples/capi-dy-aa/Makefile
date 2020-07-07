CXX = c++
CXXFLAGS = -std=c++11
PINEAPPL_DEPS != pkg-config --cflags --libs pineappl_capi
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf

dyaa: dyaa.cpp
	$(CXX) $(CXXFLAGS) $(PINEAPPL_DEPS) $(LHAPDF_DEPS) $< -o $@

PHONY: clean

clean:
	rm -f dyaa
