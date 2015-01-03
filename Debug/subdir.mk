################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../C1DGrid.cpp \
../C1DHeat.cpp \
../C2DGrid.cpp \
../C2DHeat.cpp \
../CGNUPlot.cpp \
../CGrid.cpp \
../CHeat.cpp \
../CInput.cpp \
../CLinAlg.cpp \
../CTests.cpp \
../HeatEquation.cpp 

OBJS += \
./C1DGrid.o \
./C1DHeat.o \
./C2DGrid.o \
./C2DHeat.o \
./CGNUPlot.o \
./CGrid.o \
./CHeat.o \
./CInput.o \
./CLinAlg.o \
./CTests.o \
./HeatEquation.o 

CPP_DEPS += \
./C1DGrid.d \
./C1DHeat.d \
./C2DGrid.d \
./C2DHeat.d \
./CGNUPlot.d \
./CGrid.d \
./CHeat.d \
./CInput.d \
./CLinAlg.d \
./CTests.d \
./HeatEquation.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


