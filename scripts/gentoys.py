#!/usr/bin/env python

#------------------------------------------------------------------------------
# File: likelihoodscan.py
# Created: 24-April-2013 Joe Bochenek
# $Revision$
# Description: Set of H->ZZ->4l likelihood and do a 2D likelihood likelihood
#              scan in m_4l and mu_gg 
#------------------------------------------------------------------------------


import os, sys
from ROOT import *
from time import sleep
import plot_util

# Load Poisson-Gamma RooThing
gROOT.ProcessLine(".L ../src/RooHistPoissonGamma.cxx+")


# S e t u p
# ---------------------------------------------

# Signal Strength Modifiers (define for each signal and background)
mu1 = RooRealVar("mu1","mu1",1.,0,20) 
mu2 = RooRealVar("mu2","mu2",1.,0,20) 
mu3 = RooRealVar("mu3","mu3",1.,0,20) 
mu4 = RooRealVar("mu4","mu4",1.,0,20) 
mu5 = RooRealVar("mu5","mu5",1.,0,20) 

# Input parameters
m4l = RooRealVar("m4l","m4l",115, 180) 
D = RooRealVar("bnn","bnn",0, 1) 
Dvbf = RooRealVar("bnn2","bnn2",0, 1) 

# Prepare input histograms
fin = TFile( "/home/jbochenek/work/HZZ4l_2013/scripts/hzz_templates.root")
fout = TFile( "fake_data.root", "recreate")

# Set constants
channel = "2e2mu"
era = "8TeV"
mva = "2d"

masses = ["115", "120", "123", "125", "127", "130", "135", "140", "150", "160", "170"]
stepsize = 0.2
xmin = 0.1
ymin = 0.1
channels = ["2e2mu", "4mu", "4e" ]


factor = 2.

c1 = TCanvas()

iters = 100

for mass in masses:
    print mass
    for i in range(iters):
        print "{}: {}".format(mass, i)
        for channel in channels:
            fin.cd()
            # Get qq->ZZ bkg template
            histname = "qqzz_{0}_{1}_{2}".format(channel, era, mva)
            bkg_qqzz = gDirectory.Get(histname)

            histname = "ggzz_{0}_{1}_{2}".format(channel, era, mva)
            bkg_ggzz = gDirectory.Get(histname)
            # Get the data    
            histname = "data_{0}_{1}_{2}".format(channel, era, mva)
            data = gDirectory.Get(histname)

            # ZX
            histname = "zjets_{0}_{1}".format(channel, era, mva)
            bkg_zx = gDirectory.Get(histname)    
    
            # Convert TH2s to RooDataHists and RooHistPdfs

            hbkg_qqzz =  RooDataHist("bkg_qqzz","bkg_qqzz",  RooArgList(m4l,D),   bkg_qqzz)
            rbkg_qqzz =  RooHistPdf("rbkg_qqzz", "rbkg_qqzz",   RooArgSet(m4l,D), hbkg_qqzz)

            hbkg_ggzz =  RooDataHist("bkg_ggzz","bkg_ggzz",  RooArgList(m4l,D),   bkg_ggzz)
            rbkg_ggzz =  RooHistPdf("rbkg_ggzz", "rbkg_ggzz",   RooArgSet(m4l,D), hbkg_ggzz)

            rdata  =  RooDataHist("data","data",    RooArgList(m4l,D), data)

            hbkg_zx =  RooDataHist("bkg_zx","bkg_zx",  RooArgList(m4l,D),   bkg_zx)
            rbkg_zx =  RooHistPdf("rbkg_zx", "rbkg_zx",   RooArgSet(m4l,D), hbkg_zx)

             # Get the signal
            sig_ggH = gDirectory.Get(    "ggH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )    
            sig_qqH = gDirectory.Get(    "qqH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )

            # Convert TH2s to RooDataHists and RooHistPdfs
            hsig_qqH = RooDataHist("hsig_qqH","hsig_qqH",RooArgList(m4l,D),sig_qqH)
            rsig_qqH = RooHistPdf( "rsig_qqH","rsig_qqH",RooArgSet(m4l,D) ,hsig_qqH)

            hsig_ggH = RooDataHist("hsig_ggH","hsig_ggH",RooArgList(m4l,D),sig_ggH)
            rsig_ggH = RooHistPdf( "rsig_ggH","rsig_ggH",RooArgSet(m4l,D) ,hsig_ggH)

            samples = RooArgList(rsig_qqH, rsig_ggH,  rbkg_qqzz, rbkg_ggzz, rbkg_zx)
            scales  = RooArgList(mu1,      mu2      , mu3       , mu4    , mu5  )


            #, RooArgList(mu1, mu2)
            roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(m4l, D), rdata, samples, scales, 0)

            zval = roohistpg.genHist(i)
            rdata  =  roohistpg.dataHist()

            hdata = rdata.createHistogram("fakedata_"+channel, m4l, RooFit.Binning(65), RooFit.YVar(D, RooFit.Binning(20)))
            hdata.SetName("fakedata_"+channel+"_"+mass+"_"+str(i))
    
            fout.cd()
            hdata.Write()

fout.Close()
fin.Close()