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
import histutil

histutil.setStyle()

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
D1 = RooRealVar("bnn1","bnn1",0., 1.) 
D2 = RooRealVar("bnn2","bnn2",0., 1.) 
m4l = RooRealVar("m4l","m4l",115., 180.) 

# Prepare input histograms
fin = TFile( "/home/jbochenek/work/HZZ4l_2013/scripts/hzz_templates.root")

# Set constants
channel = "2e2mu"
era = "8TeV"
mva = "cat1_2dvbf"

masses = ["122", "123", "124", "125", "126", "127", "130"]
stepsize = 0.2
xmin = 0.1
ymin = 0.1
channels = ["2e2mu", "4mu", "4e" ]
mvas = ["cat2_2dvbf", "cat1_2dvbf"]
likelihoods = []

c1 = TCanvas()



for mva in mvas:
    for channel in channels:

        likelihoods.append(TGraph2D())

        mass = "125"

        # Get qq->ZZ bkg template
        histname = "qqzz_{0}_{1}_{2}".format(channel, era, mva)
        bkg_qqzz = gDirectory.Get(histname)
        
        histname = "ggzz_{0}_{1}_{2}".format(channel, era, mva)
        bkg_ggzz = gDirectory.Get(histname)

        # Get the data    
        histname = "data_{0}_{1}_{2}".format(channel, era, mva)
        data = gDirectory.Get(histname)
    
        histname = "ggH_{0}_{1}_{2}_{3}".format(mass, channel, era, mva)
        sig_ggH = gDirectory.Get(histname)    

        histname = "qqH_{0}_{1}_{2}_{3}".format(mass, channel, era, mva)
        sig_qqH = gDirectory.Get(histname)
    
        # Convert TH2s to RooDataHists and RooHistPdfs
        hbkg_qqzz =  RooDataHist("bkg_qqzz","bkg_qqzz",  RooArgList(D1,D2),   bkg_qqzz)
        rbkg_qqzz =  RooHistPdf("rbkg_qqzz", "rbkg_qqzz",   RooArgSet(D1,D2), hbkg_qqzz)

        hbkg_ggzz =  RooDataHist("bkg_ggzz","bkg_ggzz",  RooArgList(D1,D2),   bkg_ggzz)
        rbkg_ggzz =  RooHistPdf("rbkg_ggzz", "rbkg_ggzz",   RooArgSet(D1,D2), hbkg_ggzz)

        rdata  =  RooDataHist("data","data",    RooArgList(D1,D2), data)

        hsig_qqH = RooDataHist("hsig_qqH","hsig_qqH",RooArgList(D1,D2),sig_qqH)
        rsig_qqH = RooHistPdf( "rsig_qqH","rsig_qqH",RooArgSet(D1,D2) ,hsig_qqH)

        hsig_ggH = RooDataHist( "hsig_ggH","hsig_ggH",RooArgList(D1,D2),sig_ggH)
        rsig_ggH = RooHistPdf(  "rsig_ggH","rsig_ggH",RooArgSet(D1,D2) ,hsig_ggH)

        samples = RooArgList(  rbkg_qqzz, rbkg_ggzz, rsig_qqH, rsig_ggH)
        scales  = RooArgList( mu3       , mu4, mu1,      mu2     )

        bkgfactor = 1.
        mu3.setVal(bkgfactor)
        mu4.setVal(bkgfactor)
        mu5.setVal(bkgfactor)

        null_roohistpg = RooHistPoissonGamma("null_pghistpdf","null_pghistpdf",RooArgSet(D1, D2), rdata, RooArgList(rbkg_qqzz, rbkg_ggzz), RooArgList(mu3       , mu4    , mu5  ), 0)
        nullzval = null_roohistpg.getVal()

        bkgfactor = 1.2
        mu3.setVal(bkgfactor)
        mu4.setVal(bkgfactor)
        mu5.setVal(bkgfactor)


        # Do a likelihood scan over mH and mu_gg
        for binmu_vbf in xrange(30):
            for binmu_gg in xrange(30):
                mu_gg_point  = xmin + binmu_gg * stepsize
                mu_vbf_point  = ymin + binmu_vbf * stepsize

                mu2.setVal(2*mu_gg_point)        
                mu1.setVal(2*mu_vbf_point)        

                #, RooArgList(mu1, mu2)
                roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(D1,D2), rdata, samples, scales, 0)
                zval = roohistpg.getVal() # nullzval - 
                if zval > 0:
                    likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(),  mu_gg_point,  mu_vbf_point, zval);
                else:
                    likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(),  mu_gg_point,  mu_vbf_point, zval);

        likelihoods[len(likelihoods)-1].Draw("colz")
        c1.Update()

            

#dt = plot_util.multiply_likelihood(likelihoods)
dt = plot_util.add_loglikelihood(likelihoods, 1)
dt = plot_util.setzero_likelihood(dt)

    
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
c.SaveAs("../output/nllscan_muvbf_mugg_ind.png")


    
    
c = TCanvas()
c.cd()
dt.Draw("cont2")
c.SaveAs("../output/nllscan_muvbf_mugg_contz.png")
sleep(2)

# Write stuff out
fout = TFile("../output/lh_mu_output.root", "recreate")
fout.cd()
dt.SetName("nll_mu_gg_vs_mH")
dt.SetTitle("")
dt.GetXaxis().SetTitle("#mu_{gg}")
dt.GetYaxis().SetTitle("#mu_{VBF}")
dt.GetXaxis().SetTitleOffset(0.5)
dt.GetYaxis().SetTitleOffset(0.5)
dt.GetXaxis().SetTitleSize(0.06)
dt.GetYaxis().SetTitleSize(0.06)
dt.GetHistogram().SetContour(11)
dt.Write()




c = TCanvas()
c.cd()
dt.Draw("colz")
dt.Draw("cont2 same")
c.SaveAs("../output/nllscan_muvbf_mugg.png")

fout.Close()
fin.Close()