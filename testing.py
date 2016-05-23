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
import sys

def nb_init_templates (nb_cycle, nb_amplicon_final, efficiency) :
	'''
	nb_cycle : number of cycles of the PCR
	nb_amplicon_final : number of amplicon at the end of PCR
	efficiency : calculated with PCR efficiency
	'''
	nb_init = (nb_amplicon_final) / (efficiency**nb_cycle)
	return nb_init

def actual_prog(seq, forward, reverse) :
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
	effect = round(float("".join(myListOfEfficiencies)),2)
	print str(effect)+" amplicon are synthesized from 1 template at each cycle of the PCR for the this PCR reaction"


	call(["rm", "gamResult.data"])
	call(["rm", "primerDataForR.dat"])

name = ""


if len(sys.argv) == 1 :
	print "Please input sequence, forward primer sequence, and reverse primer sequence\n ./testing.py target_sequence forward_primer_sequence reverse_sequence"

if len(sys.argv) <= 2 and len(sys.argv) > 1:
	seq = str(sys.argv[1]).upper()
	print "No primers input\nITS primers by default\nforward : TCCTCCGCTTATTGATATGC\n reverse = GGAAGTAAAAGTCGTAACAAGG"
	forward = "TCCTCCGCTTATTGATATGC"
	reverse = "GGAAGTAAAAGTCGTAACAAGG"
	actual_prog(seq,forward,reverse)

elif len(sys.argv) == 4 :
	seq = str(sys.argv[1]).upper()
	forward = str(sys.argv[2]).upper()
	reverse = str(sys.argv[3]).upper()
	actual_prog(seq,forward,reverse)

#seq = "GTAACAAGGTTTCCGTAGGTGAACATGCGGAAGGATCATTGTCGAACCCTGCAAAGCAGAACGACCCGTGAACACGTAAAAACAACCGAGCGTCGAGTGGATTTTTTGATCCACTCGATGCTCTGTCGATGCCCATTTACTTGTGTTCTTTTGGACTCGGTGGATGTGTCATTGACGCAATAACAACCCCCGGCACAATGTGTGCCAAGGAAAACTAAACTTGAGAATGCTTGTTTCATGTTGCCCCGTTCGCGGTGTGCTCTGGATTGGCTTCTTTATAATTACAAACGACTCTCGGCAACGGATATCTCGGCTCACGCATCGATGAAGAACGTAGCAAAATGCGATACTTGGTGTGAATTGCAGAATCCCGTGAACCATCGAGTTTTTGAACGCAAGTTGCGCCCGAAGCCTTTTGGCCGAGGGCACGTCTGCCTGGGCGTCACGCATCGCGTCGCCCCCAACAAATCTTTGTTGGGAGCGGATATTGGTCTCCCGTGCTCATGGCGTGGTTGGCCAAAATAGGAGTCCCTTCGATGGACGCACGAACTAGTGGTGGTCGTAAAAACCCTCGTTCTTTGTTCTGCGTTAGTCGTAAGAAAATACTCTTCAAACACCCCAATGTGTTGTCTTAGGACGACGCTTCGACCGCGACCCCAGGTCAGGCGGGACTACCCGCTGAGTTTAAGCATATCAATAA"




	#print ampliconList[0].toString()

	#foo.cssBottom()
