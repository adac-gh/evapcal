FFLAGS=  gfortran -O2
ODIR = .
all:	cascipm2 tl128 ebtable

cascnm: cascnm.o mnhlib.o cascip.o
	$(FFLAGS) -o $(ODIR)/$@ cascnm.o mnhlib.o cascip.o
cascipm2: cascipm2_128.o mnhlib.o cascip.o
	$(FFLAGS) -o $(ODIR)/$@ cascipm2_128.o mnhlib.o cascip.o
ebtable:  ebtable.o
	$(FFLAGS) -o $(ODIR)/$@  ebtable.o
tl128:          tl128.f
	$(FFLAGS) -o $(ODIR)/$@ tl128.f
tl128_mod:              tl128_mod.f
	$(FFLAGS) -o $(ODIR)/$@ tl128_mod.f
cascipm2_128.o:   cascipm2_128.f
	$(FFLAGS) -c cascipm2_128.f
mnhlib.o:         mnhlib.f
	$(FFLAGS) -c mnhlib.f
cascip.o:         cascip.f
	$(FFLAGS) -c cascip.f
ebtable.o:	ebtable.f
	$(FFLAGS) -c ebtable.f
cascnm.o:   cascnm.f
	$(FFLAGS) -c cascnm.f
clean:
	rm -f *.o