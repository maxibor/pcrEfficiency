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

from html import *
from dnaString import *
import re

class DNA(dnaString):
    """Sequence main class. Amplicon as well the Primers inherits of it."""
    
    def setLabel(self, label):
        self.label = label
        
    def setSequence(self,sequence):
        self.sequence = sequence
        self.lengthSequence = len(sequence)
    
    def getLengthSequence(self):
        """Length sequence getter."""
        
        return self.lengthSequence
        
    def getSequence(self):
        return self.sequence
    
    def cleanSequence(self):
        newLine=re.compile("%0D%0A")
        sequence2=newLine.sub('',self.sequence)

        others=re.compile("[\s\d\+]")
        self.setSequence(others.sub('',sequence2))
        
        self.degenerated()

        #handle= open("/tmp/in.pr3", "w")
        #handle.write(self.sequence)
        #handle.close()

    def degenerated(self):
        """Sequence check. Only A, C, T, G are valid letters."""
        
        alphabet=['a','c','t','g','A','C','T','G']
        setSequence=set(list(self.sequence))
        errorChar= (setSequence-set(alphabet))
        
        if len(errorChar)!=0:
            foo = html()
            foo.cssUp()
            print"""
            <div class="articles"> 
            At least one of your sequences (amplicon or primer ones) 
            has nucleotides which are diferent from 'A', 'C', 'T' and 'G'. 
            Unfortunately, this tool doesn't handle 'N' or other degenerated nucleotides. 
            The offending characteres were %s.
            <br><br>
            Please go back to the <a href="../web.html">form</a>.
            </div> 
            """ %errorChar
            
            foo.cssBottom()
            sys.exit() 
            
    def gcContent(self):
        """Calculate G+C of a given DNA string"""
        
        gcContent =float(self.sequence.count("G")+self.sequence.count("C"))/ \
            float(self.sequence.count("G")+self.sequence.count("C")+ \
            self.sequence.count("A")+ self.sequence.count("T"))
        
        self.gc = gcContent
        
        return gcContent
