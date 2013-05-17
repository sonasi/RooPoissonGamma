from ROOT import *
from time import sleep

gSystem.Load("libRooFit");
print "hi worlds"

gROOT.ProcessLine(".L ../src/RooPoissonGamma.cxx+")

# Declare observable x
x = RooRealVar("x","x",0,12) 

sigma1 = RooRealVar("sigma1","width of gaussians",1.7) 
sigma2 = RooRealVar("sigma2","width of gaussians",1.4) 
sigmalist = RooArgList(sigma1)

mean1 = RooRealVar("mean1","width of gaussians", 3.) 
mean2 = RooRealVar("mean2","width of gaussians", 2.) 
meanlist = RooArgList(mean1)

mean3 = RooRealVar("mean2","width of gaussians",3.0)

poisson = RooPoisson("pg1", "pg1", x, mean3, false)
pg2 = RooPoissonGamma("pg2", "pg2", x, meanlist, sigmalist )

tbox = TPaveText(0.5, 0.7, 0.97, 0.92, "BRNDC")
tbox.AddText("Poisson-Gamma(x|#lambda=3.0,#delta#lambda=1.7) (blue) ")
tbox.AddText("Poisson(x|#lambda=3.0) (red)")

varlist = RooArgSet(x)
x.setBins(12)

data = pg2.generate(varlist, 10000)

xframe = x.frame(RooFit.Name("xframe"),RooFit.Title(""))
data.plotOn( xframe, RooFit.LineColor(kBlack))
poisson.plotOn( xframe, RooFit.LineColor(kRed) )
pg2.plotOn( xframe, RooFit.LineColor(kBlue) )
xframe.addObject(tbox);

c1 = TCanvas()
xframe.Draw()
c1.SaveAs("poissongammaexample.png")
sleep(5)

