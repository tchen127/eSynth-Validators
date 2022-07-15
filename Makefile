#
#  This file is part of synth.
#
#  synth is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  synth is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with synth.  If not, see <http://www.gnu.org/licenses/>.
#

CC=g++ 
OPT= -Wall -O2 -fpermissive 

#Path to openbabel directory which contains the header (.h) files 
OB_INC= /usr/local/include/openbabel-2.0 

#Path to openbabel directory which contains the (.so) files 
OB_LIB= /usr/local/lib

GSL_LIB=/usr/include/gsl
GSL_INC=/usr/include/
EXE = esynth

CFLAGS= $(OPT) -I$(GSL_INC) -L$(GSL_LIB) -I$(OB_INC) -L$(OB_LIB) -lz -lgsl -lgslcblas -lopenbabel -lpthread


ODIR=./obj

DEPS = EdgeAggregator.h \
	EdgeAnnotation.h \
	Instantiator.h \
	Linker.h \
	Molecule.h \
	Brick.h \
	Utilities.h \
	Atom.h \
	ConnectableAtom.h \
	LinkerConnectableAtom.h \
	BrickConnectableAtom.h \
	Bond.h \
	AtomT.h \
	OBWriter.h \
	Options.h \
	Validator.h \
  OTFValidator.h \
  OTFValidators.h \
	obgen.h \
	Thread_Pool.h \
	IdFactory.h \
	Constants.h \
	LevelHashMap.h \
	LikeMoleculesContainer.h \
	MoleculeHashHypergraph.h \
	MinimalMolecule.h \
	SmiMinimalMolecule.h \
	FixedSortedList.h \
	EdgeDatabase.h \
	FragmentEdgeMap.h \
	SimpleFragmentGraph.h \
	zpipe.h \
	bloom_filter.hpp

_OBGEN_DEPS = obgen.h 

_WRITER_DEPS = compliantwriter.h 

_TOSMI_DEPS = convertToSMI.h 

_OBJ = Atom.o \
	Instantiator.o \
	Linker.o\
	Main.o \
	Molecule.o \
	Brick.o \
	Utilities.o \
	Atom.o \
	ConnectableAtom.o \
	LinkerConnectableAtom.o \
	BrickConnectableAtom.o \
	Bond.o \
	AtomT.o \
	OBWriter.o \
	Validator.o \
  OTFValidator.o \
  OTFValidators.o \
	Options.o \
	IdFactory.o \
	obgen.o \
	Constants.o \
	FragmentEdgeMap.o \
	zpipe.o


OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXE): $(OBJ)
	$(CC) $^ $(CFLAGS) -o $@



_OBGEN_OBJ = synthobgen.o

OBGEN_OBJ = $(patsubst %,$(ODIR)/%,$(_OBGEN_OBJ))

$(ODIR)/%.o: %.cpp $(_OBGEN_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

synthobgen: $(OBGEN_OBJ)
	$(CC) $^ $(CFLAGS) -o $@



_WRITER_OBJ = compliantwriter.o

WRITER_OBJ = $(patsubst %,$(ODIR)/%,$(_WRITER_OBJ))

$(ODIR)/%.o: %.cpp $(_WRITER_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

compliantwriter: $(WRITER_OBJ)
	$(CC) $^ $(CFLAGS) -o $@




_TOSMI_OBJ = convertToSMI.o

TOSMI_OBJ = $(patsubst %,$(ODIR)/%,$(_TOSMI_OBJ))

$(ODIR)/%.o: %.cpp $(_TOSMI_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

convertToSMI: $(TOSMI_OBJ)
	$(CC) $^ $(CFLAGS) -o $@



.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(EXE) synth.stackdump $(INCDIR)/*~
