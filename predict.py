#!/usr/bin/env python
# -*- coding: utf-8  -*-

# This file is part of PCR efficiency calculator.  
# PCR efficiency calculator is free software: you can
# redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, version 2.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Copyright Izaskun Mallona
# izaskun.mallona@gmail.com

from subprocess import call

class predict:
    """Interface to R: gam model-based prediction
    of the introduced data."""
    
    def writeHandleForR(self,myListOfAmplicons,numreturn):
        """Writes a dataframe readable by R to apply the gam model.
        The model is

            efficiency ~ s(lengthSequence, gcSequence) + s(primersLength, 

            gcPrimers) + s(gcImbalance, primerDimers).
        So a dataframe containing tab separated values of
        that variables must be generated."""
        
        handle = open('./primerDataForR.dat','w')
        i = 0
        while (i < int(numreturn)):
            amplicon = myListOfAmplicons[i]
            toWrite= str(amplicon.getLengthSequence())+ '\t' + \
                str(amplicon.gcContent()) + '\t' + \
                str(amplicon.getPrimersLength()) + '\t' + \
                str(amplicon.getGcPrimers()) + '\t' + \
                str(amplicon.getGcImbalance()) + '\t' + \
                str(amplicon.getPrimerDimers()) + ' \n'
            
            handle.write(toWrite)
            i += 1
        handle.close()

##    def efficiencySingleProduct(self):
##        """Calls to a R script which parses the data from writeHandleForR and
##        produces an output containing the estimated efficiency.
##        Only an amplicon is analyzed."""
##        
##        call('R --vanilla CMD BATCH ./gam.R', shell = True)
##        
##        handle = open('/tmp/gam.out','r')
##        efficiency = handle.read()
##        handle.close()
##        
##        return efficiency

    def predictGam(self):
        """Calls the R script and returns
        a list containing the expected efficy for each primer pair."""
        
        call('R --vanilla < gam.R > /dev/null', shell = True)
        
        handle = open('./gamResult.data','r')
        
        efficienciesList = []

        for line in handle.readlines():
            efficienciesList.append(line.rstrip())
        del efficienciesList[0]
        
        handle.close()        
        
        return efficienciesList