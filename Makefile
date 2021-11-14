
CC         = gcc
FPIC       = "-fPIC"

EXEC_NAME  =cidx

LONG_IDS   =off
PRODUCTION =on
DEBUG      =off

DEBUG_INFO =off
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
ifeq ($(DEBUG),on)
$(info internal debugging checks on)
OPT += -DDEBUG
endif

COMMON_CFLAGS = -std=c11 -fopenmp -W -Wall -I ./ $(OPT) $(FPIC)

# ----------------------------------------------------------------------------

ifeq ($(SANITIZE),on)
	SANITIZE_FLAGS = -fsanitize=float-divide-by-zero -fsanitize=bounds -fsanitize=signed-integer-overflow -fsanitize=vla-bound -fsanitize=undefined -fsanitize=integer-divide-by-zero
endif

ifeq ($(DEBUG_INFO),on)
	DBG_INFO = -ggdb3
endif

ifeq ($(PRODUCTION),on)
$(info compiling for production)
	OPTIMIZE_FLAGS = -O3 -march=native -faggressive-loop-optimizations
	CFLAGS = $(DBG_INFO) $(OPTIMIZE_FLAGS) -fno-stack-protector $(COMMON_CFLAGS) $(SANITIZE_FLAGS)
	CIDX = $(EXEC_NAME)
else
$(info compiling with debugging symbols)
	CIDX = $(EXEC_NAME)_g
	CFLAGS = -ggdb3 $(SANITIZE_FLAGS) $(COMMON_CFLAGS)
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