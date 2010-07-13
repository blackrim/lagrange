RM := rm -rf

CPP_SRCS += \
./AncSplit.cpp \
./InputReader_copper.cpp \
./BayesianBioGeoAllDispersal.cpp \
./BayesianBioGeo.cpp \
./BioGeoTree_copper.cpp \
./BioGeoTreeTools_copper.cpp \
./BranchSegment_copper.cpp \
./OptimizeBioGeo_copper.cpp \
./OptimizeBioGeoAllDispersal_copper.cpp \
./RateMatrixUtils.cpp \
./RateModel.cpp \
./Utils.cpp \
./node.cpp \
./tree.cpp \
./tree_reader.cpp \
./tree_utils.cpp \
./main_copper.cpp

OBJS += \
./AncSplit.o \
./InputReader_copper.o \
./BayesianBioGeoAllDispersal.o \
./BayesianBioGeo.o \
./BioGeoTree_copper.o \
./BioGeoTreeTools_copper.o \
./BranchSegment_copper.o \
./OptimizeBioGeo_copper.o \
./OptimizeBioGeoAllDispersal_copper.o \
./RateMatrixUtils.o \
./RateModel.o \
./Utils.o \
./node.o \
./tree.o \
./tree_reader.o \
./tree_utils.o \
./main_copper.o

CPP_DEPS += \
./AncSplit.d \
./InputReader_copper.d \
./BayesianBioGeoAllDispersal.d \
./BayesianBioGeo.d \
./BioGeoTree_copper.d \
./BioGeoTreeTools_copper.d \
./BranchSegment_copper.d \
./OptimizeBioGeo.d \
./OptimizeBioGeoAllDispersal.d \
./RateMatrixUtils.d \
./RateModel.d \
./Utils.d \
./node.d \
./tree.d \
./tree_reader.d \
./tree_utils.d \
./main_copper.d

# uncomment if debugging
# DEBUG = -DDEBUG

TARGET_NAME = lagrange_cpp

C_OPT = -O3 -ftree-vectorize -ffast-math -Weffc++ -g3
#C_OPT = -Wall -g

# for reading web input files
PYTHON_LIB = -I/usr/include/python2.6/
# output of 
# >>> import distutils.sysconfig
# >>> distutils.sysconfig.get_config_var('LINKFORSHARED')
PYTHON_REQ = -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions

#if using octave
#INCLUDES = -I/usr/include/octave-3.0.5/octave/
#INCLUDES = -I/usr/include/octave-3.2.4/

# requires fortran, gsl, and pthread -- and -ldl -lutil -lpython2.6 are for python
# if llapack lblas fail, try larmadillo
LIBS := -llapack -lblas -lgfortran -lgsl -lgslcblas -lm -lpthread -lgmp -ldl -lutil -lpython2.6

###########
# change to yes for bigtrees -- loses about 3x speed
# if 64 bit GSL try CPPFLAGS="-arch x86_64" LDFLAGS="-arch x86_64" ./configure
# need to install gmp (with ./configure --enable-cxx) and mpfr and gmpfrxx
#######
BIG = no 
BIGTREE =
ifeq  ($(strip $(BIG)),yes)
	BIGTREE += -DBIGTREE
	TARGET_NAME = lagrange_cpp_bt
	LIBS += -lgmp -lgmpxx -lmpfr -lgmpfrxx
	
endif


#######
# FORTRAN BIT
######
FC	= /usr/bin/gfortran
FFLAGS	= -O3
.f.o:;  $(FC) $(FFLAGS) -c $<

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.f
	@echo 'Building file: $<'
	@echo 'Invoking: gfortran Fortran Compiler'
	$(FC) $(FFLAGS) -c $<
	@echo 'Finished building: $<'
	@echo ' '

FORT_OBJS += \
./clock.o \
./my_expokit.o \
./mataid.o \
./blas.o \
./lapack.o \
./my_matexp.o


# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ $(DEBUG) $(BIGTREE) $(PYTHON_LIB) $(INCLUDES) $(C_OPT) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# link library locations
LINK_LIB_DIRS = -L/usr/lib/ -L/usr/local/lib/ -L./gmpfrxx/

# Tool invocations
lagrange_cpp: $(OBJS) $(FORT_OBJS)
	@echo 'Building target: $@'
	g++ $(LINK_LIB_DIRS) $(PYTHON_REQ) -o "$(TARGET_NAME)" $(FORT_OBJS) $(OBJS) $(LIBS)
	@echo 'Finished building target: $(TARGET_NAME)'
	@echo ' '

#oct: $(OBJS) $(FORT_OBJS)
#	@echo 'Building target: $@'
#	mkoctfile --link-stand-alone $(LINK_LIB_DIRS) -o "$(TARGET_NAME)" $(FORT_OBJS) $(OBJS) $(LIBS)
#	@echo 'Finished building target: $(TARGET_NAME)'
#	@echo ' '
 

# All Target
all: lagrange_cpp

# Other Targets
clean:
	-$(RM) *.o *.d $(TARGET_NAME)
	-@echo ' '
