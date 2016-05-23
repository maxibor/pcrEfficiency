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

class dnaString:
    """Common methods related to DNA string manipulation."""
    
    def setValue(self,value):
        self.value = value
        
    def complementString(self): 
        """Return the complementary sequence string.""" 

        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
        letters = list(self.value) 
        letters = [basecomplement[base] for base in letters] 
        return ''.join(letters) 

    def reverseString(self): 
        """Return the sequence string in reverse order.""" 
        
        s = self.value
        s2=s[::-1]
        return s2

    def reverseComplementString(self): 
        """Return the reverse complement of the DNA string.""" 
        
        s = self.reverseString() 
        foo = dnaString()
        foo.setValue(s)
        s2 = foo.complementString() 
        return s2