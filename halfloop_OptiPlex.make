
# Compiler
CC=gcc
CCFLAGS_ori= -fPIC -w -O2 -D_FMA -mfma -fopenmp -lm -I./dir_h
CCFLAGS_opt= -D_COMPLEX -D_AVX -mavx -D_CBLAS -Wl,-R -Wl,/usr/lib/x86_64-linux-gnu -lgslcblas
CCFLAGS_dev = $(CCFLAGS_ori) # $(CCFLAGS_opt)
CCFLAGS=$(CCFLAGS_dev) -D_MONOLITH

vpath %.c ./dir_c

project=halfloop
project_dev=halfloop_dev
project_lib=halfloop_lib
project_ar=libhalfloop

CLINK_dev=$(CC) -o $(project_dev)
CLINK=$(CC) -o $(project)

header_dev = ./dir_h/halfloop_header.h

sources_dev = array_extract_call.c \
	array_mean_center_call.c \
	array_orth_call.c \
	array_printf_call.c \
	array_stats_call.c \
	array_transpose_call.c \
	rand_call.c \
	dp_ps_call.c \
	erfcln_call.c \
	find_internal_maximum_call.c \
	get_xdrop_logscale_call.c \
	GLOBAL_tictoc.c \
	gumbel_call.c \
	halfloop_nonbinary_call.c \
	halfloop_nonbinary_recursive_call.c \
	malloc1_call.c \
	MDA_io_call.c \
	nelder_mead_call.c \
	quicksort_call.c \
	update_global_call.c \
	main_driver.c \

objects_dev = $(patsubst %.c,./dir_o/%.o,$(sources_dev)) 

all: $(objects_dev) $(project_dev)

$(objects_dev): | ./dir_o

./dir_o:
	@mkdir -p $@

$(project_dev): $(objects_dev) $(header_dev)
	rm -f $(project_dev)
	$(CLINK_dev) $(objects_dev) $(CCFLAGS_dev)

$(project): $(objects_dev) $(header_dev)
	rm -f $(project)
	$(CLINK) $(objects_dev) $(CCFLAGS)

$(project_lib): $(objects_dev) $(header_dev)
	rm -f $(project_lib).so
	gcc -shared -o $(project_lib).so $(objects_dev) $(CCFLAGS_dev)

$(project_ar): $(objects_dev) $(header_dev)
	rm -f $(project_ar).a
	ar -rcs $(project_ar).a $(objects_dev)

./dir_o/%.o : %.c ./dir_h/halfloop_header.h
	@echo $< 
	@$(CC) $(CCFLAGS_dev) -c $< -o $@

.c.o: ./dir_h/halfloop_header.h
	$(CC) -c $(CCFLAGS_dev) $<

clean: 
	rm -f $(objects_dev)

list_sources: $(sources_dev)
	echo $^

list_objects: $(objects_dev)
	echo $^

