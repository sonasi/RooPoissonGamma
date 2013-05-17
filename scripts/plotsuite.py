#!/usr/bin/env python

#------------------------------------------------------------------------------
# File: plotsuite.py
# Created: 24-April-2013 Joe Bochenek
# $Revision$
# Description: Make some plots
#------------------------------------------------------------------------------

import os, sys
from ROOT import *
from time import sleep
import plot_util
import histutil
from array import array




def flattenx(th2d):

    xPoints = lhs[0].GetX();
    yPoints = lhs[0].GetY();

    zPoints = []
    nPoints = []
    for lh in lhs:
        zPoints.append(lh.GetZ())
        nPoints.append(lh.GetN())

    dtx = []
    dty = []    
    dtz = []

    dt = TGraph2D()
    
    for j in range(nPoints[0]):
#        print j
        dtx.append(xPoints[j])
        dty.append(yPoints[j])
        zPointprod = 1
        for zPoint in zPoints:
            zPointprod *= zPoint[j]
            dt.SetPoint( j, xPoints[j], yPoints[j], zPointprod );
    return dt

histutil.setStyle()


fin1 = TFile("../output/lh_mh_output.root")
dt1 = gDirectory.Get("nll_mu_gg_vs_mH")
c1 = TCanvas()
c1.Divide(2,2)
c1.cd(1)
dt1.Draw("cont2")
c1.cd(2)
dt1.GetHistogram().ProjectionX().Draw()



ttitle=TPaveText(0.2,0.95,0.6,1.0, "bordersize=5 NDC")
ttitle.SetFillStyle(4000)
ttitle.SetFillColor(0)
ttitle.SetShadowColor(0);
ttitle.AddText("CMS Preliminary, #sqrt{s} = 8TeV, 19.6fb^{-1}")



CMSstyle.SetCanvasDefW(535) #Width of canvas
CMSstyle.SetPadRightMargin(0.17)
CMSstyle.SetPadLeftMargin(0.20)



dt1.SetTitle("")
dt1.GetXaxis().SetTitle("m_{H}")
dt1.GetYaxis().SetTitle("#mu")
dt1.GetZaxis().SetTitle("-2 #Delta ln L")
dt1.SetMaximum(20.)
c = TCanvas()
c.cd()
dt1.Draw("colz")
dt1_ = dt1.Clone()
contour = array('d', [2, 4])
dt1_.GetHistogram().SetContour(  2,  contour )
dt1_.Draw("cont2 same")
ttitle.Draw("same")
c.SaveAs("../output/nllscan_m4l_mugg.png")
c.SaveAs("../output/nllscan_m4l_mugg.pdf")



sleep(1)

fin2 = TFile("../output/lh_mu_output.root")
dt2 = gDirectory.Get("nll_mu_gg_vs_mH")
c1.cd(3)
dt2.Draw("cont2")
c1.cd(4)
dt2.GetHistogram().ProjectionX().Draw()



dt2.SetTitle("")
dt2.GetXaxis().SetTitle("#mu_{ggH}")
dt2.GetYaxis().SetTitle("#mu_{VBF}")
dt2.GetZaxis().SetTitle("-2 #Delta ln L")
dt2.SetMaximum(20.)
c = TCanvas()
c.cd()
dt2.Draw("colz")
dt2_ = dt2.Clone()
contour = array('d', [2, 4])
dt2_.GetHistogram().SetContour(  2,  contour )
dt2_.Draw("cont2 same")
ttitle.Draw("same")
c.SaveAs("../output/nllscan_mugg_muvbf.png")
c.SaveAs("../output/nllscan_mugg_muvbf.pdf")

sleep(2)


sleep(4)