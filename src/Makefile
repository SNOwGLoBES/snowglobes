# GLoBES -- General LOng Baseline Experiment Simulator
# (C) 2002 - 2007,  The GLoBES Team
#
# GLoBES is mainly intended for academic purposes. Proper
# credit must be given if you use GLoBES or parts of it. Please
# read the section 'Credit' in the README file.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# The configure script should have replaced anything as required to
# obtain a working Makefile, which may serve as a template for own
# applications.
#
# This simple Makefile is for the GLoBES examples
#
# Compile example files with ``make example1'' etc.
#
# This Makefile assumes dynamic libraries, installed to either
# the default prefix /usr/local/ or to a user-defined directory 
# called ${prefix}.
#
# For linking against a specific version of GLoBES, libglobes.so can be 
# replaced by the respective library, such as libglobes.so.0.0.1


# SNOwGLoBES compilation flags -----------------------------------------------------------------
prefix = ${GLB_DIR}
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

local_CFLAGS = -g -O4 -std=c99

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)
# ----------------------------------------------------------------------------------------------

# cevns_recoil_response.cc compilation flags ---------------------------------------------------
# These compilation flags support the Coherent Elastic Neutrino-Nucleus Scattering (CEvNS)
# smearing matrix generation code. On calling >>make install<<, libDiffSpec_1.0.a, 
# cevns_recoil_response.o, and an executable called cevns_recoil_response are built in the 
# snowglobes/smearing_code/cevns/ subdirectory. A makefile is also available in this 
# subdirectory which only compiles theCEvNS code, should you need it.

# ROOTCFLAGS and ROOTLIBS are optional arguments - some of the modules in 
# /smearing_code/dukecevns depend on ROOT. If you wish to use or develop on these modules, you 
# will need to install ROOT and compile the modules with the ROOTCFLAGS and ROOTLIBS arguments. 
# The Makefile in /smearing_code/dukecevns contains full compilation arguments for all dukecevns
# modules and allows the user to build the CEvNS code independently of SNOwGLoBES. 

# The primary subroutine in dukecevns - cevns_recoil_response - does NOT depend on ROOT; the
# basic install of SNOwGLoBES v. 1.3 will run without ROOT. ------------------------------------

#ROOTCFLAGS = `root-config --cflags`
#ROOTLIBS   = `root-config --libs`

CXXFLAGS += -I. -I./json/include -Wall -std=c++11 

LIBBASE   = DiffSpec
VER       = 1.0
TAG       = 
LIBALIAS  = $(LIBBASE)$(TAG)
LIBNAME   = $(LIBALIAS)_$(VER)

libdiffspec     = lib$(LIBNAME).a
# ----------------------------------------------------------------------------------------------


BIN = supernova
CEVNSBIN = $(libdiffspec) cevns_recoil_response
OBJ = supernova.o myio.o FormFactor.o cevns_recoil_response.o

	
%.o : %.cc
	$(RM) $@
	$(CXX) -c $(CXXFLAGS) -o $@ $*.cc

%.o : %.c
	gcc $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)


$(libdiffspec) : $(OBJ)
	$(RM) $@
	ar clq $@ $(OBJ)
	ranlib $@

formfactors: formfactors.o $(libdiffspec) 
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS) 

.PHONY: formfactors.o
formfactors.o: 
	$(CXX) -o formfactors.o $(ROOTCFLAGS) $(CXXFLAGS) -c ./../smearing_code/dukecevns/formfactors.cc

cevns_recoil_response: cevns_recoil_response.o $(libdiffspec)
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^

.PHONY: cevns_recoil_response.o $(libdiffspec)
cevns_recoil_response.o:
	$(CXX) -o cevns_recoil_response.o $(CXXFLAGS) -c ./../smearing_code/dukecevns/cevns_recoil_response.cc

supernova: supernova.o myio.o
	gcc supernova.o myio.o -o supernova  $(LDFLAGS) $(local_LDFLAGS)


all: $(BIN) $(CEVNSBIN)

install: $(BIN) $(CEVNSBIN)
	cp $(BIN) ${SNOWGLOBES}/bin
	mv $(CEVNSBIN) ${SNOWGLOBES}/smearing_code/dukecevns
	mv cevns_recoil_response.o ${SNOWGLOBES}/smearing_code/dukecevns
	mv FormFactor.o ${SNOWGLOBES}/smearing_code/dukecevns

.PHONY: clean
clean:
	rm -f $(BIN) $(CEVNSBIN) $(OBJ)

