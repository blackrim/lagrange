RM := rm -rf

CPP_SRCS += \
./AncSplit.cpp \
./InputReader.cpp \
./BayesianBioGeoAllDispersal.cpp \
./BayesianBioGeo.cpp \
./BioGeoTree.cpp \
./BioGeoTreeTools.cpp \
./BranchSegment.cpp \
./OptimizeBioGeo.cpp \
./OptimizeBioGeoAllDispersal.cpp \
./RateMatrixUtils.cpp \
./RateModel.cpp \
./Utils.cpp \
./node.cpp \
./tree.cpp \
./tree_reader.cpp \
./tree_utils.cpp \
./superdouble.cpp \
./main.cpp

OBJS += \
./AncSplit.o \
./InputReader.o \
./BayesianBioGeoAllDispersal.o \
./BayesianBioGeo.o \
./BioGeoTree.o \
./BioGeoTreeTools.o \
./BranchSegment.o \
./OptimizeBioGeo.o \
./OptimizeBioGeoAllDispersal.o \
./RateMatrixUtils.o \
./RateModel.o \
./Utils.o \
./node.o \
./tree.o \
./tree_reader.o \
./tree_utils.o \
./superdouble.o \
./main.o

CPP_DEPS += \
./AncSplit.d \
./InputReader.d \
./BayesianBioGeoAllDispersal.d \
./BayesianBioGeo.d \
./BioGeoTree.d \
./BioGeoTreeTools.d \
./BranchSegment.d \
./OptimizeBioGeo.d \
./OptimizeBioGeoAllDispersal.d \
./RateMatrixUtils.d \
./RateModel.d \
./Utils.d \
./node.d \
./tree.d \
./tree_reader.d \
./tree_utils.d \
./superdouble.d \
./main.d

# uncomment if debugging
# DEBUG = -DDEBUG

TARGET_NAME = lagrange_cpp
#for cleaning use -Weffc++
C_OPT = -O3 -ftree-vectorize -ffast-math -g3
#C_OPT = -Wall -g

# for reading web input files
#PYTHON_LIB = -I/usr/include/python2.6/
# output of 
# >>> import distutils.sysconfig
# >>> distutils.sysconfig.get_config_var('LINKFORSHARED')
#PYTHON_REQ = -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions

#if using octave
#INCLUDES = -I/usr/include/octave-3.0.5/octave/
#INCLUDES = -I/usr/include/octave-3.2.4/

# requires fortran, gsl, and pthread -- and -ldl -lutil -lpython2.6 are for python
#-ldl -lutil -lpython2.6
# if llapack lblas fail, try larmadillo
LIBS := -llapack -lblas -lgfortran -lgsl -lgslcblas -lm -lpthread -fopenmp


#######
# FORTRAN BIT
######
FC	= /usr/bin/gfortran
FFLAGS	= -O3
.f.o:;  $(FC) $(FFLAGS) -c $<

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.f
	$(FC) $(FFLAGS) -c $<

FORT_OBJS += \
./clock.o \
./my_expokit.o \
./mataid.o \
./blas.o \
./lapack.o \
./my_matexp.o


# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.cpp
	g++ $(DEBUG) $(BIGTREE) $(PYTHON_LIB) $(INCLUDES) $(C_OPT) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"

# link library locations
LINK_LIB_DIRS = -L/usr/lib/ -L/usr/local/lib/

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
	-$(RM) *.o *.d
	-@echo ' '
