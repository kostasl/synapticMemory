################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../InputVectorHandling.cpp \
../main.cpp \
../util.cpp 

OBJS += \
./InputVectorHandling.o \
./main.o \
./util.o 

CPP_DEPS += \
./InputVectorHandling.d \
./main.d \
./util.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/ -I/usr/include/boost -O3 -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


