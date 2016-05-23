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

from DNA import *
import textwrap
import re

class Amplicon(DNA):
    """An amplicon it's defined by its sequence as well the flanking primers."""
    
    def toString(self):
        """Returns the data associated to a given amplicon,
        such GC content of its primers and sequence, their
        length and the computed efficiency. """
        
        print """
        Forward: %s  (GC content: %s, length: %s) <br> 
        Reverse: %s   (GC content: %s, length: %s) <br>
        Expected efficiency: %s <br>
        Amplicon: <br>
        <tt>
        %s 
        </tt>
        <br><br>
        """ %(self.getForward() , self.forward.gcContent(), self.forward.getLengthSequence(),\
        self.getReverse(), self.reverse.gcContent(),self.reverse.getLengthSequence(),\
        self.getEfficiency(),textwrap.fill(self.getSequence(),width=60))
    
    def setEfficiency(self,efficiency):
        """Efficiency setter."""
        
        self.efficiency= efficiency

    def setPrimerPair(self,forward,reverse):
        """ Two primers conform an amplicon."""
        
        myForward = DNA()
        myForward.setLabel('forward')
        myForward.setSequence(forward)
        
        myReverse = DNA()
        myReverse.setLabel('reverse')
        myReverse.setSequence(reverse)
        
        self.forward = myForward
        self.reverse = myReverse
           
    def getGcPrimers(self):
        """Sets and returns the G+C content of both primers"""
        
        gcForward = self.forward.gcContent()
        gcReverse = self.reverse.gcContent()
        
        gcPrimers = (gcForward + gcReverse)/2
        
        self.gcPrimers = gcPrimers
        
        return self.gcPrimers
        
    
    def getGcImbalance(self):
        """ GC imbalance measure (difference between forward and reverse)"""
        
        gcForward = self.forward.gcContent()
        gcReverse = self.reverse.gcContent()
        
        gcImbalance = abs(gcForward-gcReverse)
        
        self.gcImbalance = gcImbalance
        
        return self.gcImbalance
        
    def getPrimerDimers(self):
        """Given two primer sequences a complementarity measure
        is estimated."""
        
        forwardTriplets = list()
        for x in range(len(self.forward.sequence)):
            forwardTriplets.append(self.forward.sequence[x:x+3])
        
        nudge = dnaString()
        nudge.setValue(self.reverse.sequence)
        
        reverseRevCom = nudge.reverseComplementString()
        reverseRevComTriplets = list()
        for x in range(len(reverseRevCom)):
            reverseRevComTriplets.append(reverseRevCom[x:x+3])
        
        primerDimers= int()    
        for item in forwardTriplets:
            if item in reverseRevComTriplets:
                primerDimers+=1
                
        self.primerDimers = primerDimers
        
        return primerDimers

    def getPrimersLength(self):
        """ Sums the total primer length, which is calculated
        as the sum of the individual forward and reverse
        primer sequence lengths."""
        
        forw = len(self.forward.sequence)
        reve = len(self.reverse.sequence)
        
        self.primersLength = forw + reve
        
        return self.primersLength
    
    def getForward(self):
        """Forward sequence getter."""
        
        return self.forward.sequence
    
    def getReverse(self):
        """Reverse sequence getter."""
        
        return self.reverse.sequence
    
    def getEfficiency(self):
        """Efficiency getter."""
        
        return self.efficiency

    def checkHybridization(self,sequence,forward,reverse):
        """Checks if a given pair of primers hybridizes
        within the sequence the tool is intended to calculate
        its efficiency. If false, returns a message and 
        flanks the amplicon with the forward primer and the
        complement and reversed reverse primer."""
        
        newLine=re.compile("%0D%0A")
        forward2=newLine.sub('',forward)

        others=re.compile("[\s\d\+]")
        forward=others.sub('',forward2)
        
        reverse2=newLine.sub('',reverse)
        reverse=others.sub('',reverse2)

        #if primers are not included in the sequence, they should be
        check=0
        forwd = re.compile(forward,re.IGNORECASE)
        revse=re.compile(reverse,re.IGNORECASE)
        found=0

        newFor = dnaString()
        newRev = dnaString()
        newSeq = dnaString()
        
        newFor.setValue(forward)
        newRev.setValue(reverse)
        newSeq.setValue(sequence)
        
        forward = newFor
        reverse = newRev
        sequence = newSeq
        
        if forwd.search(sequence.value) is not None:
            if revse.search(sequence.reverseComplementString()) is not None:
                start=sequence.value.find(forward.value)
                end=sequence.value.find(reverse.reverseComplementString())
                found=1
                
        else:
            if forwd.search(sequence.reverseComplementString()) is not None:
                if revse.search(sequence.value) is not None:
                    start=sequence.value.find(forward.value)
                    end=sequence.value.find(reverse.reverseComplementString())
                    found=1
            
        if found:
            sequence.value= sequence.value[start:end+len(reverse.value)]
            
        else:
            check=1
            sequence.value = forward.value+sequence.value
            #this is a logical element to advice the absence of primer annealing into the amplicon
            sequence.value = sequence.value + reverse.reverseComplementString() 
            
        return check, sequence.value
            
            