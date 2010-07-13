RM := rm -rf

CPP_SRCS += \
./AncSplit.cpp \
./InputReader.cpp \
./BioGeoTree.cpp \
./BioGeoTreeTools.cpp \
./BranchSegment.cpp \
./OptimizeBioGeo.cpp \
./OptimizeBioGeoAllDispersal.cpp \
./OptimizeBioGeoPowell.cpp \
./RateMatrixUtils.cpp \
./RateModel.cpp \
./Utils.cpp \
./main.cpp

OBJS += \
./AncSplit.o \
./InputReader.o \
./BioGeoTree.o \
./BioGeoTreeTools.o \
./BranchSegment.o \
./OptimizeBioGeo.o \
./OptimizeBioGeoAllDispersal.o \
./OptimizeBioGeoPowell.o \
./RateMatrixUtils.o \
./RateModel.o \
./Utils.o \
./main.o

CPP_DEPS += \
./AncSplit.d \
./InputReader.d \
./BioGeoTree.d \
./BioGeoTreeTools.d \
./BranchSegment.d \
./OptimizeBioGeo.d \
./OptimizeBioGeoAllDispersal.d \
./OptimizeBioGeoPowell.d \
./RateMatrixUtils.d \
./RateModel.d \
./Utils.d \
./main.d

C_OPT = -arch i386 -O3 -ftree-vectorize -ffast-math -Wall -g
#C_OPT = -arch i386 -O3 -Wall -g


# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Users/smitty/programming/LIBRARIES/boost_1_39_0 -I/System/Library/Frameworks/Python.framework/Headers/ $(C_OPT) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

#requires fortran, bio++ libraries, gsl, and pthread
LIBS := -lgfortran -lbpp-popgen -lbpp-phyl -lbpp-seq -lbpp-numcalc -lbpp-utils -lgsl -lgslcblas -lm -lpthread

FC	= /usr/local/bin/gfortran
FFLAGS	= -O3
.f.o:;  $(FC) $(FFLAGS) -c $<

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.f
	@echo 'Building file: $<'
	@echo 'Invoking: g77 Fortran Compiler'
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


# All Target
all: lagrange_cpp

# Tool invocations
lagrange_cpp: $(OBJS) $(FORT_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X C++ Linker'
	g++ -L/usr/local/lib/ -L/usr/lib/ -u _PyMac_Error /System/Library/Frameworks/Python.framework/Versions/2.6/Python -arch i386 -o "lagrange_cpp" $(FORT_OBJS) $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) *.o *.d lagrange_cpp
	-@echo ' '
