
#Uncomment this line to generate an optimized executable
#STD    = -c -r8 -i4 -O1
#STD    = -c -r8 -i4 -fastsse -Mnolre
#Uncomment this line to generate a debug executable
#STD    =  -c -i4  -tp k8-64 -g -Mbounds -Ktrap=fp 
#STD    =  -c -i4 -tp k8-64 -g -Mbounds 
COMP   = ftn -g -c $(STD)
#Uncomment this line to generate an optimized executable
#COMP1  = ftn -c -r8 -i4 -fastsse -Mnolre
#COMP1  =  ftn -c -r8 -i4 -O1
#Uncomment this line to generate a debug executable
#COMP1  = ftn -c -i4 -tp k8-64 -g -Mbounds -Ktrap=fp
#COMP1  = ftn -c -i4 -tp k8-64 -g -Mbounds 
RadHyd_LIB =$(HDF) $(LAPACK)
#RadHyd_LIB =$(HDF)
SEGLDR = ftn -g -o hdf5Benchmark
EXE = hdf5Benchmark
default : $(EXE)



OBJFILES=   KIND_MODULE.o  \
            NUMERICAL_MODULE.o  \
            data_module.o \
            IO_MODULE.o \
            hdf5Benchmark.o

%.o: $(VPATH)%.f95
	$(COMP) $< ${HDF5_INC}

%.o: $(VPATH)%.f
	$(COMP1) $< 

hdf5Benchmark: $(OBJFILES)
#	$(COMP) $(STD) $(OBJFILES)
	$(SEGLDR) $(RadHyd_LIB) $(OBJFILES) 
#	cp hdf5Benchmark /tmp/work/reubendb/workspace/hdf5Benchmark/
#	rm -f /tmp/work/reubendb/workspace/hdf5Benchmark/*.h5

debug:
	@$(MAKE) WHAT=_RadHyd_GR COMP="$(COMP1) $(DBG)"

clean:
	rm -f *.o *.mod
clobber: clean
	rm -f hdf5Benchmark
