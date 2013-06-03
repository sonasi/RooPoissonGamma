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
import glob

histutil.setStyle()

from array import array

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

# Set constants
channel = "2e2mu"
era = "8TeV"
mva = "2d"

masses = ["118", "120", "122", "123", "124", "126", "127", "128", "130", "135"]
massint = map(int, masses)
stepsize = 0.2
xmin = 0.1
ymin = 0.1
points = 25
channels = [ "4mu", "4e", "2e2mu" ]

likelihoods1d = []

# Prepare input histograms
findata = TFile( "fake_data.root")
fin = TFile( "/home/jbochenek/work/HZZ4l_2013/scripts/hzz_templates.root")
fout = TFile("../output/lh_mh_output_toys_2.root", "recreate")

factor = 2.


nsyst = 2

likelihoods_iter = []

itoys = 25

for iter in range(itoys):
    print "Iteration: {0}/{1}".format( iter, itoys )

    likelihoods = []

    for channel in channels:
        # Prepare input histograms
        mass = 125

        fin.cd()
        likelihoods.append(TGraph2D())

        # Get qq->ZZ bkg template

        histname = "qqzz_{0}_{1}_{2}".format(channel, era, mva)
        print histname
        bkg_qqzz = gDirectory.Get(histname)
        
        histname = "ggzz_{0}_{1}_{2}".format(channel, era, mva)
        print histname
        bkg_ggzz = gDirectory.Get(histname)

        # ZX
        histname = "zjets_{0}_{1}".format(channel, era, mva)
        bkg_zx = gDirectory.Get(histname)
        fin.cd()

        # Get the data    
        findata.cd()
        histname = "fakedata_{0}_{1}_1".format(channel, mass, iter)
        print histname
        data = gDirectory.Get(histname)
        fin.cd()


        
        

        print "mass: {0}, plot: {1}, channel: {2}".format(mass, mva, channel)

        sig_ggH = gDirectory.Get( "ggH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )    
        sig_qqH = gDirectory.Get( "qqH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )


    
        # Convert TH2s to RooDataHists and RooHistPdfs

        hbkg_qqzz =  RooDataHist("bkg_qqzz","bkg_qqzz",  RooArgList(m4l,D),   bkg_qqzz)
        rbkg_qqzz =  RooHistPdf("rbkg_qqzz", "rbkg_qqzz",   RooArgSet(m4l,D), hbkg_qqzz)

        hbkg_ggzz =  RooDataHist("bkg_ggzz","bkg_ggzz",  RooArgList(m4l,D),   bkg_ggzz)
        rbkg_ggzz =  RooHistPdf("rbkg_ggzz", "rbkg_ggzz",   RooArgSet(m4l,D), hbkg_ggzz)


        hbkg_zx =  RooDataHist("bkg_zx","bkg_zx",  RooArgList(m4l,D),   bkg_zx)
        rbkg_zx =  RooHistPdf("rbkg_zx", "rbkg_zx",   RooArgSet(m4l,D), hbkg_zx)


        rdata  =  RooDataHist("data","data",    RooArgList(m4l,D), data)

        null_roohistpg = RooHistPoissonGamma("null_pghistpdf","null_pghistpdf",RooArgSet(m4l, D), rdata, RooArgList(rbkg_qqzz, rbkg_ggzz, rbkg_zx), RooArgList(mu3       , mu4    , mu5  ), 0)
        nullzval = null_roohistpg.getVal()


        # Get the signal
        sig_ggH = gDirectory.Get(    "ggH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )    
        sig_qqH = gDirectory.Get(    "qqH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )
        hsig_qqH = RooDataHist("hsig_qqH","hsig_qqH",RooArgList(m4l,D),sig_qqH)
        rsig_qqH = RooHistPdf( "rsig_qqH","rsig_qqH",RooArgSet(m4l,D) ,hsig_qqH)

        hsig_ggH = RooDataHist("hsig_ggH","hsig_ggH",RooArgList(m4l,D),sig_ggH)
        rsig_ggH = RooHistPdf( "rsig_ggH","rsig_ggH",RooArgSet(m4l,D) ,hsig_ggH)

        samples = RooArgList(rsig_qqH, rsig_ggH,  rbkg_qqzz, rbkg_ggzz, rbkg_zx)
        scales  = RooArgList(mu1,      mu2      , mu3       , mu4    , mu5  )

        roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(m4l, D), rdata, samples, scales, 0)

        zval = roohistpg.genHist(iter)
        rdata  =  roohistpg.dataHist()



    #    nulllh = RooHistPoissonGamma("nullpghistpdf","nullpghistpdf",RooArgSet(m4l, D), rdata, RooArgList(rbkg_qqzz, rbkg_ggzz), RooArgList(mu1), 0)

        scanmu1d = [0]*len(masses)

        for i, mass in enumerate(masses):
                print "Channel: {0}, Mass: {1}".format(channel, mass)
    
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

                bkgfactor = 1.
                mu3.setVal(bkgfactor)
                mu4.setVal(bkgfactor)
                mu5.setVal(bkgfactor)

                # Do a likelihood scan over mH and mu_gg
                for binmu_gg in xrange(points):
                    mu_gg_point  = xmin + binmu_gg * stepsize
                    mu2.setVal(mu_gg_point)        
                    mu1.setVal(mu_gg_point)        

                    #, RooArgList(mu1, mu2)
                    roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(m4l, D), rdata, samples, scales, 0)
                    zval = 2*roohistpg.getVal() # 2*(nullzval - roohistpg.getVal())
                    if zval > 0:
                        likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(), float(mass), mu_gg_point, zval)

                mu2.setVal(2)        
                mu1.setVal(2)        
                zval = 2*roohistpg.getVal() # 2*(nullzval - roohistpg.getVal())
                scanmu1d[i] += zval

        for i, mass in enumerate(masses):
            scanmu1d[i] = scanmu1d[i]/points
        print len(masses)
        print len(massint)
        print massint
        print scanmu1d    
        gr = TGraph(len(masses), array("d", massint) , array("d",scanmu1d) )
        gr.SetName("nll_mu_mH_{0}".format(channel))
        fout.cd()
        gr.Write()
        likelihoods1d.append(gr)
    
    tot_lh = plot_util.add_loglikelihood(likelihoods, 1)
    tot_lh.SetName("nll_mu_gg_vs_mH_{0}".format(iter))
    fout.cd()
    tot_lh.Write()
    likelihoods_iter.append(plot_util.add_loglikelihood(likelihoods, 1))

#dt = plot_util.multiply_likelihood(likelihoods)
dt = plot_util.ave_ln_loglikelihood(likelihoods_iter)
#reset to zero 
dt = plot_util.setzero_likelihood(dt)

c1 = TCanvas()
c1.Divide(2)

c1.cd(1)
likelihoods_iter[0].Draw("colz")
c1.cd(2)
dt.Draw("colz")
c1.Update()

dt1d = plot_util.add_tgraph(likelihoods1d)
dt1d.SetName("nll_mH_1d")
dt1d.Write()

dt1d.Draw()




# Write stuff out
fout.cd()
dt.SetName("nll_mu_gg_vs_mH")
dt.SetTitle("")
dt.GetXaxis().SetTitle("m_{4l}")
dt.GetYaxis().SetTitle("#mu_{gg}")
#dt.GetXaxis().SetTitleOffset(0.5)
#dt.GetYaxis().SetTitleOffset(0.5)
#dt.GetXaxis().SetTitleSize(0.06)
#dt.GetYaxis().SetTitleSize(0.06)
dt.GetHistogram().SetContour(  25  )
dt.Write()




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
dt.Draw("colz")
dt.Draw("cont2 same")
c.SaveAs("../output/nllscan_m4l_mugg_contz.png")




c = TCanvas()
c.cd()
dt.Draw("colz")
c.SaveAs("../output/nllscan_m4l_mugg.png")




fout.Close()
fin.Close()



