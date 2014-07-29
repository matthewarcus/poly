APPS := 
APPS11 := poly twod schwarz

all: $(APPS) $(APPS11)

OPT := -O1

$(APPS): %: %.cpp
	g++ -MMD $(OPT) -g -Wall -o $@ $< $(EXTRA)

$(APPS11): %: %.cpp
	g++ --std=c++11 -MMD $(OPT) -g -Wall -Wshadow -o $@ $< $(EXTRA)

poly poly2 animate tut2 tut4: EXTRA += -lGL -lGLU -lglut

test: lnew2
	./lnew2 -1 1 1 0.5 0.01 0.005 0.0002 | ./animate

test2: lnew2
	./lnew2 -1 1 1 0 0 0.02 0.001 | ./animate

clean:
	rm -f $(APPS) $(APPS11) *.d

.PHONY: clean test

-include *.d
