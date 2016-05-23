#!/usr/bin/env R
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

# loading the library
suppressMessages(library(mgcv))

# loading the model
load('gam.RData')

dat<-read.table('./primerDataForR.dat',sep='\t',header=F)

# this is based on
#
# miniGaussian<-(gam(efficiency ~ s(lengthSequence,gcSequence)+s(primersLength,gcPrimers)+s(gcImbalance,primerDimers), data=dataMini))
#
#    a = query.getLengthSequence()
#    b = query.gcContent()
#    c = query.getPrimersLength()
#    d = query.getGcPrimers() 
#    e = query.getGcImbalance()
#    f = query.getPrimerDimers()

colnames(dat)<-c('lengthSequence','gcSequence','primersLength','gcPrimers', 'gcImbalance','primerDimers')

efficiencies <- data.frame()

efficiencies <- predict(miniGaussian,dat)

efficiencies[as.numeric(efficiencies)<1.5]<-'1.5'
efficiencies[as.numeric(efficiencies)>2]<-'2'


write.table(efficiencies, './gamResult.data', row.names = FALSE,quote=FALSE)
