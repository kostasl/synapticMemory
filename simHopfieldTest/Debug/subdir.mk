################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../HopfieldMemoryTests.cpp \
../InputVectorHandling.cpp \
../main.cpp \
../synapseAllocators.cpp \
../util.cpp 

OBJS += \
./HopfieldMemoryTests.o \
./InputVectorHandling.o \
./main.o \
./synapseAllocators.o \
./util.o 

CPP_DEPS += \
./HopfieldMemoryTests.d \
./InputVectorHandling.d \
./main.d \
./synapseAllocators.d \
./util.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/home/kostasl/CodeProjects/synapticMemory/synapseModels" -I/usr/include/boost -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


