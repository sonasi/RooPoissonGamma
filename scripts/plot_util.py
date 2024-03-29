from ROOT import *
import array, os, sys, re

def cdf(hist):
	n = hist.GetNbinsX()
	c = []
	print n
	sum = 0
	for i in range(n):
		sum = sum + hist.GetBinContent(i)
		c.append(sum)
	return c;
	


#------------------------------------------------------------------------------
# Multiply some tgraphs together
#------------------------------------------------------------------------------

def multiply_likelihood(lhs):

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


#------------------------------------------------------------------------------
# Multiply some tgraphs together
#------------------------------------------------------------------------------

def add_loglikelihood(lhs, i):

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
            zPointprod += zPoint[j]
            dt.SetPoint( j, xPoints[j], yPoints[j], float(zPointprod)/float(i) );
    return dt



#------------------------------------------------------------------------------
# Set the minimum of the likelihood to zero ( so subract the min from all values)
#------------------------------------------------------------------------------

def setzero_likelihood(lh):

    xPoints = lh.GetX();
    yPoints = lh.GetY();
    zPoints = lh.GetZ();
    
    nPoints = lh.GetN()
    zmin = lh.GetZmin()
    dt = TGraph2D()
    print "Likelihood Minimum: {}".format(zmin)
    
    for j in range(nPoints):
        print "{}/{} {}".format(j, nPoints, zPoints[j])
        zPointnew = zPoints[j] - zmin
        dt.SetPoint( j, xPoints[j], yPoints[j], zPointnew );
    return dt



#------------------------------------------------------------------------------
# Multiply some tgraphs together
#------------------------------------------------------------------------------

def add_tgraph(lhs):

    xPoints = lhs[0].GetX();
    zPoints = []

    nPoints = []
    for lh in lhs:
        zPoints.append(lh.GetY())
        nPoints.append(lh.GetN())

    dtx = []
    dty = []    

    dt = TGraph()
    
    for j in range(nPoints[0]):
#        print j
        dtx.append(xPoints[j])
        zPointprod = 1
        for zPoint in zPoints:
            zPointprod += zPoint[j]
            dt.SetPoint( j, xPoints[j], float(zPointprod) );
    return dt



#------------------------------------------------------------------------------
# Calculate an estimate of the significance using a likelihood and null LH
#------------------------------------------------------------------------------
def significance(lh, nulllh):

    xPoints = lh.GetX();
    yPoints = lh.GetY();
    zPoints = lh.GetZ()
    nPoints = lh.GetN()

    dtx = []
    dty = []    
    dtz = []

    dt = TGraph2D()
    
    for j in range(nPoints):
        zPointprod = zPoints[j]/nulllh
        signif = 0
        if zPoints[j] > nulllh:
            signif = sqrt(2*log(zPoints[j]/nulllh)) 
        dt.SetPoint( j, xPoints[j], yPoints[j], signif );
    return dt

#------------------------------------------------------------------------------
# Unfold Contents 
#------------------------------------------------------------------------------
def unfoldContents(h):
	c = []
	for i in range(h.GetNbinsX()):
		binx = i+1
		for j in range(h.GetNbinsY()):
			biny = j+1
			c.append( h.GetBinContent(binx, biny) )
	return c

#------------------------------------------------------------------------------
# Make Discriminant distributions (loop thorugh events)
#------------------------------------------------------------------------------

def loop(ifile, hbnn, mva, cmva, mvavars, type):
    print ifile
    filein = TFile(ifile)
    treein = filein.Get("HZZ4LeptonsAnalysis")
    nentry = treein.GetEntries()

    count = 0 	

	#do signal	
    for event in treein:
        count+=1

        weight = float(event.f_weight)
        if event.f_njets_pass < 2:
            continue
        inputs = []
        for var in mvavars:
            inputs.append(getattr(treein, var))

        if(type==1):
            y = mva(inputs, 50, 100)
            hbnn.Fill(y)        
        if(type==0):
            y = mva(inputs)
            hbnn.Fill(y)

        if( not (count%1000)):
            print count
            cmva.cd()
#            hbnn.Draw()
            cmva.Update()
        
    filein.Close()	
    print "Total Count: " + str(count)	



#------------------------------------------------------------------------------
# Make individual ROC plot
#------------------------------------------------------------------------------

def calceff(geff, heff, hs, hb):

    stotal=hs.Integral();
    btotal=hb.Integral();

    cs = cdf(hs);
    cb = cdf(hb);

    maxsig = 0;
    optmvax = 0;
    optmvay = 0;
    optcut = 0;

    bins = hs.GetNbinsX()

    for i in range(hs.GetNbinsX()):
        cs1 = stotal - cs[i-1]
        cb1 = btotal - cb[i-1]
        es = cs1 / stotal
        eb = cb1 / btotal

        print "bins:{}, es: {}, eb: {} ".format(i, es, eb)
        sign=0
        if(cb1 > 0): 
            sign = cs1/sqrt(cs1 + cb1)
        if(sign > maxsig):
            maxsig = sign;
            optcut = i;	 
            optmvax = eb;
            optmvay = es;
            heff.Fill(eb, es);
        if(i<2):
            continue
        else:
            geff.SetPoint(geff.GetN(), eb, es);

        

    print "\n Maximum Significance " + str(maxsig) + ", is bin: " + str(optcut) + ", which is BNN Cut: " + str(optcut/bins)
    print "es: " + str(optmvay) + "\teb: " + str(optmvax)





def getrate(file, process):
    in_file = open(file)
    for line in in_file:
        data = line.rsplit()
        if data:
            if data[0] == "rate" and data[1] == process:
                return data[2]
    return -1

def getefficiency(file, ratefile, dataset):
    in_file = open(ratefile)
    count = getcount(file)
    for line in in_file:
        data = line.rsplit()
        if data:
            if data[0] == dataset:
                return float(count)/float(data[1])
    return -1

def getcs(file, ratefile, dataset):
    in_file = open(ratefile)
    for line in in_file:
        data = line.rsplit()
        if data:
            if data[0] == dataset:
                return float(data[4])
    return -1
    
def getcount(file):
    in_file = open(file)
    line = in_file.readline()
    data = line.rsplit()
    columns = []
    for varname in data:
        varname = varname.lstrip(':b')
        columns.append(varname); 
    count = 0

    for line in in_file:
        data = line.rsplit()
        d = dict(zip(columns, data))
        count += 1
    return count-1
    
def getsumweight(files):
    weight_sum = 0
    for file in files:
        in_file = open(file)
        line = in_file.readline()
        data = line.rsplit()
        columns = []
        for varname in data:
            varname = varname.lstrip(':b')
            columns.append(varname)
        for line in in_file:
            data = line.rsplit()
            d = dict(zip(columns, data))
            weight_sum += float(d["weight"])
    return weight_sum
    
def getsumweight_ntuple(files):
    weight_sum = 0
    for file in files:
        print file
        f = TFile(file)
        tree = f.Get("HZZ4LeptonsAnalysis")
        # Fill Histograms with tight and loose candidates
        for event in tree:
            if 100 < event.f_mass4l < 1000:
                weight_sum += event.f_weight
    return weight_sum
