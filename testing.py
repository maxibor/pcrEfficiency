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


#this avoids zero float division errors

from __future__ import division
from subprocess import call
from amplicon import *
from predict import *
import operator


name = "test"
seq = "GTAACAAGGTTTCCGTAGGTGAACATGCGGAAGGATCATTGTCGAACCCTGCAAAGCAGAACGACCCGTGAACACGTAAAAACAACCGAGCGTCGAGTGGATTTTTTGATCCACTCGATGCTCTGTCGATGCCCATTTACTTGTGTTCTTTTGGACTCGGTGGATGTGTCATTGACGCAATAACAACCCCCGGCACAATGTGTGCCAAGGAAAACTAAACTTGAGAATGCTTGTTTCATGTTGCCCCGTTCGCGGTGTGCTCTGGATTGGCTTCTTTATAATTACAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGCCGAGGGCACGTCTGCCTGGGCGTCACGCATCGCGTCGCCCCCAACAAATCTTTGTTGGGAGCGGATATTGGTCTCCCGTGCTCATGGCGTGGTTGGCCAAAATAGGAGTCCCTTCGATGGACGCACGAACTAGTGGTGGTCGTAAAAACCCTCGTTCTTTGTTCTGCGTTAGTCGTAAGAAAATACTCTTCAAACACCCCAATGTGTTGTCTTAGGACGACGCTTCGACCGCGACCCCAGGTCAGGCGGGACTACCCGCTGAGTTTAAGCATATCAATAA"
reverse = "GGAAGTAAAAGTCGTAACAAGG"
forward = "TCCTCCGCTTATTGATATGC"


query = Amplicon()
query.setLabel(name)
query.setSequence(seq)

#check for non [ACTG] chars
query.cleanSequence()

forward = forward.upper()
reverse = reverse.upper()

query.setSequence(query.sequence.upper())
query.setPrimerPair(forward,reverse)

#check for non [ACTG] chars
query.cleanSequence()

#foo = html()
#foo.cssUp()

checking = query.checkHybridization(query.getSequence(),query.forward.sequence,query.reverse.sequence)
query.setSequence(checking[1])

	
ampliconList = []
ampliconList.append(query)

#query.checkHybridization()
myPredict = predict()
myPredict.writeHandleForR(ampliconList,1)
myListOfEfficiencies = myPredict.predictGam()

ampliconList[0].setEfficiency(myListOfEfficiencies[0])
#print name
effect = int("".join(myListOfEfficiencies))
print str(effect)+" amplicon are synthesized from 1 template at each cycle of the PCR for the this '"+name+"' PCR reaction"


call(["rm", "gamResult.data"])
call(["rm", "primerDataForR.dat"])


#print ampliconList[0].toString()

#foo.cssBottom()


