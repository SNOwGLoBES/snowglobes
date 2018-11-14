# Let's find a parameterization which does a good job of matching "known" cross sections
# This code works on cross section shape
# code initially written and tested with python 2.7.10 and ROOT 6.12
# It may not work for python 3.x due to the ROOT import
# by jba 8/9/2018

import ROOT
import code
import array
import scipy.optimize as optimize
import math

# no shadowboxes!
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)


# first, a function to read in a cross section file and create a histogram
def ReadXS(filename, flavornum = 1, xsname = "unnamed"):
   # 1 = nue, 2 = numu, 3 = nutau, 4 = nuebar, etc.
   xs_data = []
   infile = open(filename, "r")
   for line in infile.readlines():
      if line[0] == "#":
         continue # skip over header lines
      if len(line) < 5:
         continue # skip over empty lines
      else:
         split_line = line.split()
         xs_data.append([float(x) for x in split_line])
   # now we make a histogram
   numbins = len(xs_data)
   logstep = (xs_data[30][0] - xs_data[10][0])/20.0
   minbinlowedge = (xs_data[0][0] - logstep/2.0)
   binlist = []
   for ibin in range (numbins+1):
      binlist.append(minbinlowedge + ibin * logstep)
   binarray = array.array("d", [pow(10, x) for x in binlist])
   maxbinhighedge = (xs_data[-1][0] + logstep/2.0)
   flavornamelist = ["nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]
   hist = ROOT.TH1D("xs_" + xsname + "_" + flavornamelist[flavornum-1], 
                     xsname + "_" + flavornamelist[flavornum-1] + ";E (GeV);10^{-38} cm^{2}",
                     numbins, binarray)
   for xs_datum in xs_data:
      hist.Fill(pow(10, xs_datum[0]), xs_datum[flavornum]*pow(10, xs_datum[0]))
   for ibin in range(numbins):
      hist.SetBinError(ibin+1, 0)
   return hist

c12hist = ReadXS("../../snowglobes/xscns/xs_nue_C12.dat", 1, "C12cc")
ar40hist = ReadXS("../../snowglobes/xscns/xs_nue_Ar40.dat", 1, "Ar40cc")
dhist = ReadXS("../../snowglobes/xscns/xs_nue_d.dat", 1, "dcc")
c13hist = ReadXS("../../snowglobes/xscns/xs_nue_C13.dat", 1, "C13cc")
o16hist = ReadXS("../../snowglobes/xscns/xs_nue_O16.dat", 1, "O16cc")

# Now, we need to try to make a function which fits all the cross sections well

def FitFunc(x, q, a, p0, p1, p2, p3):
   xp = x-q
   if x < q:
      return 0
   val = 1000*xp*xp*xp + (p0 + p1/a + q*1000*p2)*xp*xp + p3/a*xp*xp*xp*xp
   # below are some functional forms that were tried and rejected
   #val = 1000*xp*xp + (p0 + math.sqrt(a)*p1)*xp + p3
   #val = xp + p1 * xp*xp + p2*a*pow(xp, 2 + p3) + (p4 + p5*math.sqrt(a))*pow(xp, 2 + p0)
   #val = xp + p1 * xp*xp + (p2 + p3*math.sqrt(a)) * x*x + (p4 + p5*math.sqrt(a))*pow(xp, p0)
   #val = (p5 + math.sqrt(a))*xp*xp*xp + (p1+p2*math.sqrt(a)) * xp + (p3+p4*math.sqrt(a)) * xp*xp
   #val = ((p0*a)*pow(xp, 1+p4*a) + p1*xp*xp + p2*(pow(xp,1+p3*a)))
   return val


# A chi-square sort of thing, but ONLY FOR SHAPE
# The histograms are normalized to each other before testing
def GoodnessOfFit(q, a, p0, p1, p2, p3, hist):
   nbins = hist.GetNbinsX()
   chi2thingy = 0
   nbinsfit = 0
   histsum = hist.Integral()
   fitsum  = 0
   for ibin in range(nbins):
      binx = hist.GetBinCenter(ibin+1)
      fitsum += FitFunc(binx, q, a, p0, p1, p2, p3)
   for ibin in range(nbins):
      binx = hist.GetBinCenter(ibin+1)
      biny = hist.GetBinContent(ibin+1)
      fity = FitFunc(binx, q, a, p0, p1, p2, p3)*histsum/fitsum
      if biny > 0:
         chi2thingy += pow(fity-biny,2)/biny * 10000 + math.sqrt(pow(fity-biny,2))/biny * 100
         nbinsfit += 1
      else:
         chi2thingy += 0
   chi2thingy /= nbinsfit
   return chi2thingy

# the total goodness of fit over all the histos
def TotalGoodness(params):
   p0, p1, p2, p3 = params
   histlist = [[c12hist, 18.84*0.001, 12],
               [o16hist, 19.76*0.001, 16],
               [ar40hist, 4.3*0.001, 40],
               [dhist, 1.44*0.001, 2]]
   totalchi2thingy = 0
   for histpair in histlist:
      totalchi2thingy += GoodnessOfFit(histpair[1], histpair[2], p0, p1, p2, p3, histpair[0])
   print (params, totalchi2thingy)
   return totalchi2thingy


# Ok, now we can just do a simple optimization!
initial_guess = [26.25, 153.35,  -1.934,  1000]
bounds_list = [(-5000, 5000), (-5000, 5000), (-5000, 5000), (-5000, 50000)]
result = optimize.minimize(TotalGoodness, initial_guess, method = "SLSQP",
                           bounds = bounds_list,
                           options = {"maxiter":10000})
print(result)

def FitFuncHist(sample_hist, name, params, q, a):
   p0, p1, p2, p3 = params
   newhist = sample_hist.Clone("h_" + name + "_cloned")
   newhist.Reset()
   nbins = newhist.GetNbinsX()
   for ibin in range(nbins):
      binx = newhist.GetBinCenter(ibin+1)
      biny = FitFunc(binx, q, a, p0, p1, p2, p3)
      newhist.Fill(binx, biny)
   for ibin in range(nbins):
      newhist.SetBinError(ibin, 0)
   return newhist

histlist = [[c12hist, 18.84*0.001, 12],
            [o16hist, 19.76*0.001, 16],
            [ar40hist, 4.3*0.001, 40],
            [dhist, 1.44*0.001, 2]]



canvaslist = [0]*4
fittedhistlist = [0]*4

for i, histpair in enumerate(histlist):
   canvaslist[i] = ROOT.TCanvas()
   ihist = histpair[0]
   ihist.Draw()
   fittedhistlist[i] = FitFuncHist(ihist, ihist.GetName(), result.x, histpair[1], histpair[2])
   fittedhistlist[i].SetLineColor(ROOT.kMagenta)
   fittedhistlist[i].Scale(ihist.Integral()/fittedhistlist[i].Integral())
   fittedhistlist[i].Draw("same")
#   canvaslist[i].SetLogx()
#   canvaslist[i].SetLogy()



code.interact("let's go!", local = globals())




