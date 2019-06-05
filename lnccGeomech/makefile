# the compiler to use
compiler	:=	gfortran	

# compiler flags:
cflags	= -c -Wall -Wconversion -g -fbacktrace -fbounds-check -ffpe-trap=zero,overflow -ffree-line-length-512 -Og
lflags	= -g

# program name
program	= plasprog

# list of directories
srcFolder	=./fonte
binFolder		=./bin
objFolder		=./bin

# source files and objects
objsfiles	= variaveisGlobais.o malha.o leituraEscrita.o funcoesDeForma.o utilitarios.o mInputReader.o propGeoFisica.o solverGaussSkyline.o leituraEscritaSimHidroGeoMec.o geomecanica.o hidrodinamica.o driverGeocreep.o
fullObjs	= $(patsubst %,$(objFolder)/%,$(objsfiles))

all: $(binFolder)/$(program)

$(binFolder)/$(program): $(fullObjs)
	$(compiler) $(lflags) $(fullObjs) -o $(binFolder)/$(program)

$(objFolder)/variaveisGlobais.o:

$(objFolder)/malha.o: \
	$(objFolder)/variaveisGlobais.o

$(objFolder)/leituraEscrita.o: \
	$(objFolder)/malha.o

$(objFolder)/funcoesDeForma.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/leituraEscrita.o

$(objFolder)/utilitarios.o: \
	$(objFolder)/variaveisGlobais.o

$(objFolder)/mInputReader.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/malha.o

$(objFolder)/propGeoFisica.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/malha.o

$(objFolder)/solverGaussSkyline.o:

$(objFolder)/leituraEscritaSimHidroGeoMec.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/malha.o \
	$(objFolder)/leituraEscrita.o \
	
$(objFolder)/geomecanica.o: \
	$(objFolder)/variaveisGlobais.o \

$(objFolder)/hidrodinamica.o:

$(objFolder)/driverGeocreep.o:



$(objFolder)/%.o: $(srcFolder)/%.F90
	$(compiler) $(cflags) -o $@ $<

