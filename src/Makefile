include ../Makefile.inc

SHARED = libPoisonGamma.so
INCLUDELIST= RooPoissonGamma.h

all: RooPoissonGamma rootdict.cc $(SHARED)

rootdict.cc: $(INCLUDELIST)
	@echo $(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) $^ 
	@echo $(CXXFLAGS)
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) $^ 

RooPoissonGamma: 
	@echo $(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c RooPoissonGamma.cxx 
	@echo "CXX Flags "
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c RooPoissonGamma.cxx 

$(SHARED):  RooPoissonGamma.o 
	@echo "---> Linking $@"
	@echo $(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^
	mv $(SHARED) ../lib

clean:
	@rm -f *.o
	@rm -f rootdict.*
	@rm -f $(SHARED)
	@rm -f ../lib/$(SHARED)
