
CC = g++
CFLAGS = -Wall -c -g
LFLAGS = -Wall -g

comp: main.o
	$(CC)  $(LFLAGS)  -o  comp Direct.o  Alignment.o  Matrix.o  wrapper.o  main.o  Encoded.o
# Each command must begin with a TAB.

main.o: main.cc  Alignment.o  Direct.o  Encoded.o  Matrix.o  wrapper.o 
	$(CC)  $(CFLAGS)  main.cc

Alignment.o: Alignment.cc  Alignment.hh  Direct.hh Matrix.hh wrapper.hh
	$(CC)  $(CFLAGS)  Alignment.cc

Direct.o: Direct.cc  Direct.hh wrapper.hh
	$(CC)  $(CFLAGS)  Direct.cc

Encoded.o: Encoded.cc Encoded.hh Alignment.hh wrapper.hh Direct.hh
	$(CC)  $(CFLAGS) Encoded.cc

wrapper.o: wrapper.cc wrapper.hh
	$(CC)  $(CFLAGS)  wrapper.cc

#template.stack.min.o: template.stack.min.hh
#	$(CC)  $(CFLAGS)  template.stack.min.hh

Matrix.o: Matrix.cc Matrix.hh Alignment.hh Direct.hh wrapper.hh
	$(CC)  $(CFLAGS)  Matrix.cc



# Removes all .o files and executables.
# The backslash '\' stops "rm" from prompting the user
# for whether to remove the file.

#clean : 
#	\rm  *.o  comp
