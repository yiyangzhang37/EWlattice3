### This makefile simply includes the shortcuts for mpirun.
# detect OS
ifeq ($(OS), Windows_NT)
	detected_OS := Windows
else
	detected_OS := $(shell uname -s)
endif

MPICXX = mpic++
MPIRUN = mpirun

TARGET = emlat.out
OUTPUT = result.txt

run:$(TARGET)
	$(MPIRUN) -np 1 ./$(TARGET)
run1:$(TARGET)
	$(MPIRUN) -np 1 ./$(TARGET) -r 1 -c 1 &> $(OUTPUT) &
run2:$(TARGET)
	$(MPIRUN) -np 2 ./$(TARGET) -r 1 -c 2 &> $(OUTPUT) &
run4:$(TARGET)
	$(MPIRUN) -np 4 ./$(TARGET) -r 2 -c 2 &> $(OUTPUT) &
run8:$(TARGET)
	$(MPIRUN) -np 8 ./$(TARGET) -r 2 -c 4 &> $(OUTPUT) &
run32:$(TARGET)
	$(MPIRUN) -np 32 ./$(TARGET) -r 4 -c 8 &> $(OUTPUT) &
run32:$(TARGET)
	$(MPIRUN) -np 40 ./$(TARGET) -r 5 -c 8 &> $(OUTPUT) &

clean:
	$(RM) $(OBJS) $(TARGET)
cleandata:
	$(RM) *.h5 result.txt *_dtable.txt *_param.txt