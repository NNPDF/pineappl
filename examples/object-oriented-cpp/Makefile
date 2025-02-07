CXX = c++
CXXFLAGS = -std=c++11 -g -Wall -Wextra -O3
PINEAPPL_DEPS != pkg-config --cflags --libs pineappl_capi
LHAPDF_DEPS != pkg-config --cflags --libs lhapdf

PROGRAMS = \
	dyaa \

all: $(PROGRAMS)

test-examples: $(PROGRAMS)
	set -e && for i in $(PROGRAMS); do ./$${i}; done > test-output && diff -u output test-output && rm test-output

dyaa: dyaa.cpp
	$(CXX) $(CXXFLAGS) $< $(PINEAPPL_DEPS) $(LHAPDF_DEPS) -o $@

PHONY: clean

clean:
	rm -f $(PROGRAMS) *.pineappl *.lz4
