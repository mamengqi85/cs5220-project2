include Makefile.in

.PHONY: all exe doc clean realclean run

exe: sph.x 
doc: main.pdf derivation.pdf
all: exe doc

# =======

sph.x: sph.o params.o state.o particle.o bin.o interact.o leapfrog.o io_bin.o dump.o timing.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

sph.o: sph.c params.h state.h particle.o bin.o interact.h leapfrog.h io.h dump.h timing.h

params.o: params.c params.h
state.o: state.c state.h
particle.o: particle.c particle.h dump.h
bin.o: bin.c bin.h particle.h dump.h
interact.o: interact.c interact.h state.h particle.h bin.h params.h dump.h
leapfrog.o: leapfrog.c leapfrog.h state.h particle.h bin.h params.h dump.h
dump.o: dump.c dump.h state.h particle.h bin.h
io_txt.o: io_txt.c io.h particle.h bin.h
io_bin.o: io_bin.c io.h particle.h bin.h

%.o: %.c
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $<

io_bin.o: io_bin.c
	$(CC) -c $(CFLAGS) $<

# =======
main.pdf: main.tex codes.tex
derivation.pdf: derivation.tex check_derivation.tex

codes.tex: params.h state.h particle.h interact.c leapfrog.c sph.c params.c io_bin.c
	dsbweb -o $@ -c $^

check_derivation.tex: check_derivation.m
	dsbweb -o $@ -m $^

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

# =======
view: 
	java -jar ../jbouncy/Bouncy.jar run.out

# =======
tgz: realclean
	(cd ..; tar -czf sph.tgz sph)

clean:
	rm -f *~ *.o
	rm -f codes.tex check_derivation.tex
	rm -f main.log main.aux main.out main.toc
	rm -f derivation.log derivation.aux derivation.out

realclean: clean
	rm -f main.pdf derivation.pdf *.x run.out

# ======
run: sph.x
	qsub run.qsub

viewer:
	java -jar Bouncy.jar run.out

diff:
	diff run.out ../sph_origin/run.out
