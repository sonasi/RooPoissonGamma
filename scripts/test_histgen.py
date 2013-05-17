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

# Parameters
mu1 = RooRealVar("mu1","mu1",1.,0,10) 
mu2 = RooRealVar("mu2","mu2",1.,0,10) 

m4l = RooRealVar("m4l","m4l",120) 
D = RooRealVar("bnn","bnn",0.5) 

# Prepare input histograms
fin = TFile("hzz_templates.root")
fin.cd()

# Set constants
channel = "2e2mu"
era = "8TeV"
mva = "cat2_2d"

masses = ["122", "123", "124", "125", "126", "127", "130"]
stepsize = 0.9
xmin = 0.2
ymin = 0.1
channels = ["4e", "4mu", "2e2mu"]

likelihoods = []


for channel in channels:
    likelihoods.append(TGraph2D())

    # Get qq->ZZ bkg template
    histname = "zz_train_" + channel + "_" + era + "_" + mva
    zzbkg = gDirectory.Get(histname)

    # Get the data
    histname = "data_{0}_{1}_{2}".format(channel, era, mva)
    data = gDirectory.Get(histname)

    # Convert TH2s to RooDataHists and RooHistPdfs
    rdata  =  RooDataHist("data","data",    RooArgList(m4l,D),data)
    hzzbkg =  RooDataHist("zzbkg","zzbkg",  RooArgList(m4l,D),zzbkg)
    rzzbkg =  RooHistPdf("rzzbkg", "rzzbkg",   RooArgSet(m4l,D), hzzbkg)

    nulllh = RooHistPoissonGamma("nullpghistpdf","nullpghistpdf",RooArgSet(m4l, D), rdata, RooArgList(rzzbkg), RooArgList(mu1), 0)

    mass = "123"
    sig = gDirectory.Get(    "sig_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )    
    hsignal = RooDataHist("signal","signal",RooArgList(m4l,D),sig)
    rsignal = RooHistPdf("rsignal", "rsignal", RooArgSet(m4l,D), hsignal)
    roohistpg1 = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(m4l, D), rdata, RooArgList(rzzbkg), RooArgList(mu1), 0)

    zval = roohistpg1.genHist(1223)
    rdata  =  roohistpg1.dataHist()
    hzzbkg =  RooDataHist("zzbkg","zzbkg",  RooArgList(m4l,D),zzbkg)
    rzzbkg =  RooHistPdf("rzzbkg", "rzzbkg",   RooArgSet(m4l,D), hzzbkg)

    print roohistpg1.dataHist().createHistogram("this", m4l, RooFit.Binning(20), RooFit.YVar(D, RooFit.Binning(20))).Draw("colz")
    sleep(1)
    
    print zzbkg.Integral()

    for mass in masses:
            print "Channel: {}, Mass: {}".format(channel, mass)
    
            # Get the signal
            sig = gDirectory.Get(    "sig_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )    
            sigvbf = gDirectory.Get(  "sig_VBF_" + str(mass) + "_" + channel + "_"+era+"_" + mva    )

            # Convert TH2s to RooDataHists and RooHistPdfs
            hsignal = RooDataHist("signal","signal",RooArgList(m4l,D),sig)
            rsignal = RooHistPdf("rsignal", "rsignal", RooArgSet(m4l,D), hsignal)

            # Do a likelihood scan over mH and mu_gg
            for binmu_gg in xrange(20):
                mu_vbf_point = 1.
                mu_gg_point  = xmin + binmu_gg * stepsize

                mu2.setVal(mu_vbf_point)
                mu1.setVal(mu_gg_point)        

                #, RooArgList(mu1, mu2)
                mu2.setVal(1.0)
                roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(m4l, D), rdata, RooArgList(rzzbkg), RooArgList(mu1), 0)
                zval = roohistpg.getVal()
                likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(), float(mass), mu_gg_point, zval);

    likelihoods[len(likelihoods)-1].Draw("colz")
    sleep(1)
            

#dt = plot_util.multiply_likelihood(likelihoods)
dt = plot_util.add_loglikelihood(likelihoods)

c = TCanvas()
c.Divide(2,2)
c.cd(1)
likelihoods[0].Draw("colz")
c.cd(2)
likelihoods[1].Draw("colz")
c.cd(3)
likelihoods[2].Draw("colz")
c.cd(4)
dt.Draw("colz")

c.SaveAs("../output/nllscan_m4l_mugg_ind.png")


    
c = TCanvas()
c.cd()
dt.Draw("cont2")
c.SaveAs("../output/nllscan_m4l_mugg_contz.png")
sleep(2)

# Write stuff out
fout = TFile("../output/lh_output.root", "recreate")
fout.cd()

dt.SetName("nll_mu_gg_vs_mH")
dt.SetTitle("")
dt.GetXaxis().SetTitle("m_{4l}")
dt.GetYaxis().SetTitle("#mu_{gg}")

dt.Write()

c = TCanvas()
c.cd()
dt.Draw("colz")
c.SaveAs("../output/nllscan_m4l_mugg.png")



fout.Close()
fin.Close()