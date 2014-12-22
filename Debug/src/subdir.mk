################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/list.c \
../src/parallel_selection.c 

OBJS += \
./src/list.o \
./src/parallel_selection.o 

C_DEPS += \
./src/list.d \
./src/parallel_selection.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	mpicc -I/usr/lib/mpich/include -I/usr/local/bin -I/usr/include/mpich -I"/root/workspace/parallel_selection/src" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


