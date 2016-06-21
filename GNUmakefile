THISDIR:=$(abspath $(dir $(realpath $(lastword ${MAKEFILE_LIST}))))
CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
LDLIBS+= -lm -lpthread
LDFLAGS+= $(patsubst %,-L%,${LIBDIRS})
ALL = datander TANmask REPmask HPC.TANmask HPC.REPmask
vpath %.c ${THISDIR}

%: %.c

all: ${ALL}
datander: tandem.o align.o
TANmask: align.o
REPmask: align.o
${ALL}: align.o DB.o QV.o

clean:
	rm -f ${ALL}
	rm -fr *.dSYM *.o

install:
	rsync -av ${ALL} ${PREFIX}/bin
symlink:
	ln -sf $(addprefix ${CURDIR}/,${ALL}) ${PREFIX}/bin

.PHONY: clean all
