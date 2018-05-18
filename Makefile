DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = datander TANmask REPmask HPC.TANmask HPC.REPmask HPC.DAScover

all: $(ALL)

datander: datander.c tandem.c tandem.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o datander datander.c tandem.c align.c DB.c QV.c -lpthread -lm

TANmask: TANmask.c align.h align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o TANmask TANmask.c align.c DB.c QV.c -lm

REPmask: REPmask.c align.h align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o REPmask REPmask.c align.c DB.c QV.c -lm

HPC.TANmask: HPC.TANmask.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPC.TANmask HPC.TANmask.c DB.c QV.c -lm

HPC.REPmask: HPC.REPmask.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPC.REPmask HPC.REPmask.c DB.c QV.c -lm

HPC.DAScover: HPC.DAScover.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPC.DAScover HPC.DAScover.c DB.c QV.c -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f damasker.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf damasker.tar.gz README.md Makefile *.h *.c
