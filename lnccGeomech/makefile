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

# source files and create a list of all object files
objsfiles	= variaveisGlobais.o malha.o leituraEscrita.o funcoesDeForma.o utilitarios.o mInputReader.o propGeoFisica.o solverGaussSkyline.o leituraEscritaSimHidroGeoMec.o geomecanica.o hidrodinamica.o driverGeocreep.o
fullObjs	= $(patsubst %,$(objFolder)/%,$(objsfiles))

# make without parameters default to make all
all: $(binFolder)/$(program)

# compile for all objects
$(objFolder)/%.o: $(srcFolder)/%.F90
	$(compiler) $(cflags) -J$(objFolder) -o $@ $<

# Linker
$(binFolder)/$(program): $(fullObjs)
	$(compiler) $(lflags) $(fullObjs) -o $(binFolder)/$(program)

# Debug command to show the variables
debug:
	@echo "objsfiles = $(objsfiles)"
	@echo "fullObjs = $(fullObjs)"

clean:
	rm -f $(objFolder)/*.o $(objFolder)/*.mod $(binFolder)/$(program)

.PHONY: debug all clean

# Dependencies
$(objFolder)/driverGeocreep.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/malha.o \
	$(objFolder)/propGeoFisica.o \
	$(objFolder)/leituraEscrita.o \
	$(objFolder)/leituraEscritaSimHidroGeoMec.o \
	$(objFolder)/utilitarios.o \
	$(objFolder)/geomecanica.o \
	$(objFolder)/hidrodinamica.o

$(objFolder)/geomecanica.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/funcoesDeForma.o \
	$(objFolder)/malha.o \
	$(objFolder)/propGeoFisica.o \
	$(objFolder)/leituraEscritaSimHidroGeoMec.o \
	$(objFolder)/solverGaussSkyline.o \
	$(objFolder)/utilitarios.o

$(objFolder)/hidrodinamica.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/funcoesDeForma.o \
	$(objFolder)/malha.o \
	$(objFolder)/propGeoFisica.o \
	$(objFolder)/solverGaussSkyline.o \
	$(objFolder)/leituraEscritaSimHidroGeoMec.o \
	$(objFolder)/utilitarios.o

$(objFolder)/solverGaussSkyline.o: \
	$(objFolder)/utilitarios.o

$(objFolder)/leituraEscritaSimHidroGeoMec.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/malha.o \
	$(objFolder)/mInputReader.o \
	$(objFolder)/leituraEscrita.o \
	$(objFolder)/propGeoFisica.o

$(objFolder)/funcoesDeForma.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/utilitarios.o \
	$(objFolder)/leituraEscrita.o

$(objFolder)/malha.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/utilitarios.o

$(objFolder)/leituraEscrita.o: \
	$(objFolder)/malha.o

$(objFolder)/mInputReader.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/malha.o

$(objFolder)/propGeoFisica.o: \
	$(objFolder)/variaveisGlobais.o \
	$(objFolder)/malha.o

$(objFolder)/utilitarios.o: \
	$(objFolder)/variaveisGlobais.o

$(objFolder)/variaveisGlobais.o: