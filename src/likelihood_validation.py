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

# Prepare input histograms
fin = TFile( "/home/jbochenek/work/HZZ4l_2013/scripts/hzz_templates.root")
fin2 = TFile( "fake_data.root")

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

iters = 25

likelihoods = []

for i in range(iters):
    print "ITER: " + str(i)
    for channel in channels:
        likelihoods.append(TGraph2D())

        fin.cd()
        # Get qq->ZZ bkg template
        histname = "qqzz_{0}_{1}_{2}".format(channel, era, mva)
        bkg_qqzz = gDirectory.Get(histname)

        histname = "ggzz_{0}_{1}_{2}".format(channel, era, mva)
        bkg_ggzz = gDirectory.Get(histname)

        fin2.cd()
        # Get the data    
        histname = "fakedata_{0}_{1}".format(channel, i, era, mva)
        data = gDirectory.Get(histname)
        fin.cd()

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



        mass = "126"
        sig = gDirectory.Get(    "ggH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )


        stack = THStack()
        stack.Add(bkg_zx.ProjectionX())
        stack.Add(bkg_ggzz.ProjectionX())
        stack.Add(bkg_qqzz.ProjectionX())

        stack.Add(sig.ProjectionX())    
        sig.ProjectionX().SetFillColor(kRed)
        sig.ProjectionX().SetFillStyle(3001)
        sig.ProjectionX().SetLineColor(kRed)

        bkg_qqzz.ProjectionX().SetFillColor(kGreen)
        bkg_qqzz.ProjectionX().SetLineColor(kGreen)

        bkg_qqzz.ProjectionX().SetFillColor(kYellow)
        bkg_qqzz.ProjectionX().SetLineColor(kYellow)

        bkg_zx.ProjectionX().SetFillColor(kBlue)
        bkg_zx.ProjectionX().SetLineColor(kBlue)

        null_roohistpg = RooHistPoissonGamma("null_pghistpdf","null_pghistpdf",RooArgSet(m4l, D), rdata, RooArgList(rbkg_qqzz, rbkg_ggzz, rbkg_zx), RooArgList(mu3       , mu4    , mu5  ), 0)
        nullzval = null_roohistpg.getVal()


    #    nulllh = RooHistPoissonGamma("nullpghistpdf","nullpghistpdf",RooArgSet(m4l, D), rdata, RooArgList(rbkg_qqzz, rbkg_ggzz), RooArgList(mu1), 0)

        for mass in masses:
                print "Channel: {}, Mass: {}".format(channel, mass)
    
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

                # Do a likelihood scan over mH and mu_gg
                for binmu_gg in xrange(20):
                    mu_gg_point  = xmin + binmu_gg * stepsize
                    mu2.setVal(mu_gg_point)        

                    #, RooArgList(mu1, mu2)
                    roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(m4l, D), rdata, samples, scales, 0)

                    zval = nullzval - roohistpg.getVal() 
                    if zval > 0:
                        likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(), float(mass), mu_gg_point, zval);
                    else:
                        likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(), float(mass), mu_gg_point, zval);

    #dt = plot_util.multiply_likelihood(likelihoods)
    dt = plot_util.add_loglikelihood(likelihoods, int(i+1))

    print int(dt.GetZmax() - dt.GetZmin())
    dt.GetHistogram().SetContour(int(dt.GetZmax() - dt.GetZmin()))

    dt.Draw("colz")   
    dt.Draw("cont2 same")

    c1.Update()         

    
    
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

c.SaveAs("../output/nllscan_m4l_mugg_fake_ind.png")


    
    
c = TCanvas()
c.cd()
dt.Draw("cont2")
c.SaveAs("../output/nllscan_m4l_mugg_fake_contz.png")
sleep(2)

# Write stuff out
fout = TFile("../output/lh_output_fake.root", "recreate")
fout.cd()

dt.SetName("nll_mu_gg_vs_mH")
dt.SetTitle("")
dt.GetXaxis().SetTitle("m_{4l}")
dt.GetYaxis().SetTitle("#mu_{gg}")
dt.GetXaxis().SetTitleOffset(0.5)
dt.GetYaxis().SetTitleOffset(0.5)
dt.GetXaxis().SetTitleSize(0.06)
dt.GetYaxis().SetTitleSize(0.06)
dt.GetHistogram().SetContour(int(dt.GetZmax() - dt.GetZmin()))
dt.Write()

c = TCanvas()
c.cd()
dt.Draw("colz")
dt.Draw("cont2 same")
c.SaveAs("../output/nllscan_m4l_mugg_fake.png")

fout.Close()
fin.Close()