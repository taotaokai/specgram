############################################################################
# 'A Generic Makefile for Building Multiple main() Targets in $PWD'
# Author:  Robert A. Nader (2012)
# Email: naderra at some g
# Web: xiberix
############################################################################
#  The purpose of this makefile is to compile to executable all C source
#  files in CWD, where each .c file has a main() function, and each object
#  links with a common LDFLAG.
#
#  This makefile should suffice for simple projects that require building
#  similar executable targets.  For example, if your CWD build requires
#  exclusively this pattern:
#
#  cc -c $(CFLAGS) main_01.c
#  cc main_01.o $(LDFLAGS) -o main_01
#
#  cc -c $(CFLAGS) main_2..c
#  cc main_02.o $(LDFLAGS) -o main_02
#
#  etc, ... a common case when compiling the programs of some chapter,
#  then you may be interested in using this makefile.
#
#  What YOU do:
#
#  Set PRG_SUFFIX_FLAG below to either 0 or 1 to enable or disable
#  the generation of a .exe suffix on executables
#
#  Set CFLAGS and LDFLAGS according to your needs.
#
#  What this makefile does automagically:
#
#  Sets SRC to a list of *.c files in PWD using wildcard.
#  Sets PRGS BINS and OBJS using pattern substitution.
#  Compiles each individual .c to .o object file.
#  Links each individual .o to its corresponding executable.
#
###########################################################################
#
I := ../inc
CFLAGS_INC := -I$I
#CFLAGS := -g -Wall -Wextra -pedantic $(CFLAGS_INC)
CFLAGS := -O2 -Wall -Wextra -pedantic $(CFLAGS_INC)

LDFLAGS := -lpar -lsac -lm -lfftw3 -lnetcdf

CHEADER := $I/func
PROG := ../bin/specgram
#
## ==================- NOTHING TO CHANGE BELOW THIS LINE ===================
##
SRCS := $(wildcard *.c)
## OBJS are automagically compiled by make.
OBJS := $(patsubst %.c,%.o, $(SRCS))
##
#all: $(CHEADER) $(PROG)
all: $(PROG)
##
##
#$(CHEADER): $(SRCS)
#	cfunctions --backup -n -g $(CHEADER) $(SRCS)
##
## Compile the executables
$(PROG) : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(PROG)
##
## $(OBJS) should be automatically removed right after linking.
##
clean:
	$(RM) $(PROG) $(OBJS)
##
rebuild: clean all
##
## eof Generic_Multi_Main_PWD.makefile
