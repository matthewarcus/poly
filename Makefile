APPS := poly

all: $(APPS)

OPT := -O1

$(APPS): %: %.cpp
	g++ --std=c++11 -MMD $(OPT) -g -Wall -Wshadow -o $@ $< $(EXTRA)

poly: EXTRA += -lGL -lGLU -lglut

clean:
	rm -f $(APPS) *.o *.d

.PHONY: clean

-include *.d
