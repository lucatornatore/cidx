#myCC=gcc
#myCC=clx

arch := $(shell uname -p)
$(info ${arch} architecture found)
ifeq ($(arch),ppc64le)
	myCC=xlc
else 
	myCC=gcc
endif

CC=$(myCC)
$(info using $(CC) c compiler)
FPIC       = "-fPIC"

EXEC_NAME  =cidx

LONG_IDS   =off
PRODUCTION =on

DEBUG_INFO =on
SANITIZE   =off
PAPI       =off

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

OPT=

#OPT+= -DUSE_LIBC_QSORT      # if you want to use the stadc version of qsort (slower)
#OPT+= -DUSE_LIBC_BSEARCH    # if you want to use the stadc version of bsearch (slower)

ifeq ($(LONG_IDS),on)
$(info using 8bytes ids)
OPT += -DLONG_IDS
else
$(info using 4bytes ids)
endif
ifeq ($(PAPI),on)
$(info using PAPI)
OPT += -DPAPI_MEASURE         # to obtain some metrics
endif
ifeq ($(PRODUCTION),off)
DEBUG=on
endif


ifeq ($(CC),gcc)
	OPENMP = -fopenmp 
else ifeq ($(CC),xlc)
	OPENMP = -qsmp=omp
endif

COMMON_CFLAGS = -std=c11 $(OPENMP) -W -Wall -I ./ $(OPT) $(FPIC)

# ----------------------------------------------------------------------------

ifeq ($(SANITIZE),on)
	SANITIZE_FLAGS = -fsanitize=float-divide-by-zero -fsanitize=bounds -fsanitize=signed-integer-overflow -fsanitize=vla-bound -fsanitize=undefined -fsanitize=integer-divide-by-zero
endif

ifeq ($(DEBUG_INFO),on)
	DBG_INFO = -ggdb3
endif

ifeq ($(PRODUCTION),on)
$(info compiling for production)
	ifeq (($CC),gcc)
	OPTIMIZE_FLAGS = -O4 -march=native -faggressive-loop-optimizations
	else ifeq ($(CC),xlc)
	OPTIMIZE_FLAGS = -O5 -qnostrict -qarch=auto -qtune=auto -qhot
	endif
	CFLAGS = $(DBG_INFO) $(OPTIMIZE_FLAGS) -fno-stack-protector $(COMMON_CFLAGS) $(SANITIZE_FLAGS)
	CIDX = $(EXEC_NAME)
else
$(info internal debugging checks on)
	OPT += -DDEBUG
$(info compiling with debugging symbols)
	CIDX = $(EXEC_NAME)_g
	CFLAGS = -O2 -ggdb3 $(SANITIZE_FLAGS) $(COMMON_CFLAGS)
endif

# ------------------------------------ define dependence files for each target
OBJS   = cidx_vars.o cidx.o io.o sort_and_search.o subfind_utils.o
INCL   = cidx.h mybsearch.h qsort_template.h


# ------------------------------------ define rules

cidx: $(OBJS) Makefile
	@$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(CIDX)
	@echo "========================  ::  "$(CIDX)" built"
	@echo "                              with flags " $(CFLAGS)
	@echo
	@cat *.o.warnings > warnings
	@rm -f *.o.warnings
	@echo "pls check whether there are messages in the 'warnings' file"


$(OBJS): $(INCL) Makefile

%.o : %.c
	@echo "........................[$(CC)]" $< "->" $@
	@$(CC) -c $(CFLAGS) $< -o $@ 2> $@.warnings


clean:
	@echo "cleaning objects files: "	
	@echo $(OBJS) $(CIDX) 
	@rm -f $(OBJS) $(CIDX) warnings

