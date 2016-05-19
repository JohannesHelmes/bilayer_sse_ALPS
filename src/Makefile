#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#
#  Makefile for ALPS code on local thp-machine
#       Johannes Helmes


ALPS_ROOT = /scratch/alps

include $(ALPS_ROOT)/share/alps/include.mk

CXXFLAGS = -O2 -DNDEBUG -march=corei7

all: XYbilayer


XYbilayer.o: XYbilayer.cpp XYbilayer.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++0x -c XYbilayer.cpp -o XYbilayer.o $(COMFLAG)
XYbilayer: XYbilayer.o main.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ main.o XYbilayer.o  $(LDFLAGS) $(LIBS) $(LINKFLAG)

main.o: main.C
	$(CXX) $(CPPFLAGS)  $(CXXFLAGS) -std=c++0x -c main.C


clean:
	$(RM) $(RMFLAGS) *.o XYbilayer
