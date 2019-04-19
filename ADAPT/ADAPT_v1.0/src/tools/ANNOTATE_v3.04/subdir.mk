################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ANNOTATE_v3.04.cpp \
../annotate_sequence.cpp \
../matrix_properties.cpp \
../read_matrix.cpp \
../region_read.cpp \
../score_calc.cpp 

OBJS += \
./ANNOTATE_v3.04.o \
./annotate_sequence.o \
./matrix_properties.o \
./read_matrix.o \
./region_read.o \
./score_calc.o 

CPP_DEPS += \
./ANNOTATE_v3.04.d \
./annotate_sequence.d \
./matrix_properties.d \
./read_matrix.d \
./region_read.d \
./score_calc.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


