# Quick, plot a couple cross sections against each other!
# Will work with python 2.6, but possibly not 3.x

import ROOT
import array
import sys
import math
import code

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)

if len(sys.argv) < 3:
   print ("Usage: python ./compare_xs.py <filename1> <filename2> [<filename3>...] [flavornum]\n")
   print ("flavornum is 1 = nue, 2 = numu, 3 = nutau, 4 = nuebar, 5 = numubar, 6 = nutaubar\n")
   sys.exit()


flavornum = 1
filelist = []
for ifile in range(len(sys.argv)-1):
   print sys.argv[ifile+1]
   if str(sys.argv[ifile+1]).isdigit():
      flavornum = int(sys.argv[ifile+1])
      continue
   else:
      filelist.append(sys.argv[ifile+1])

# quick code to skip yellow, the worst color in ROOT
def SkipYellow(number):
   if number > 4:
      return number + 1
   else:
      return number

# first, a script to read in a cross section file and create a histogram
def ReadXS(filename, flavornum = 1, xsname = "unnamed"):
   # 1 = nue, 2 = numu, 3 = nutau, 4 = nuebar, etc.
   print flavornum
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

c1 = ROOT.TCanvas()

histlist = []
for i, ifile in enumerate(filelist):
   histlist.append(ReadXS(ifile, flavornum, "xs_" + str(i)))

for i, ihist in enumerate(histlist):
   ihist.SetLineColor(SkipYellow(i+1))
   if i == 0:
      ihist.Draw()
   else:
      ihist.Draw("same")

code.interact("Let's go!", local=globals())
