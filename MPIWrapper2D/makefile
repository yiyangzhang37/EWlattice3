MPICXX = mpic++
MPIRUN = mpirun

# define any compile-time flags
CFLAGS = -std=c++11

# define any directories containing header files other than /usr/include
#
INCLUDES = 

# define the C source files
MPI_SRCS = MPIWrapper2D.cpp
TEST_SRCS = Test_MPIWrapper2D.cpp
SRCS = main.cpp $(MPI_SRCS) $(TEST_SRCS)

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.cpp=.o)

# define the executable file 
TARGET = lat_test
OUTPUT = lat_test.txt

all:$(TARGET)
	@echo Compile success.

$(TARGET): $(OBJS)
	$(MPICXX) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS)
.cpp.o:
	$(MPICXX) $(CFLAGS) $(INCLUDES) -c $< -o $@
run:$(TARGET)
	$(MPIRUN) -np 6 ./$(TARGET) -r 3 -c 2 > $(OUTPUT) &
run1:$(TARGET)
	$(MPIRUN) -np 1 ./$(TARGET) -r 1 -c 1 > $(OUTPUT) &
run4:$(TARGET)
	$(MPIRUN) -np 4 ./$(TARGET) -r 2 -c 2 > $(OUTPUT) &
run8:$(TARGET)
	$(MPIRUN) -np 8 ./$(TARGET) -r 2 -c 4 > $(OUTPUT) &
run15:$(TARGET)
	$(MPIRUN) -np 15 ./$(TARGET) -r 5 -c 3 > $(OUTPUT) &
clean:
	$(RM) $(OBJS) $(TARGET)