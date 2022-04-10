#!/usr/bin/env python3

import sys,os,time,getopt,requests,urllib,re,math,operator
import subprocess,pickle
import json

from itertools import combinations
from itertools import combinations_with_replacement
from itertools import permutations

### Need to make sure these libraries can be loaded! ###
import math,random
import xyplotter,webbrowser
import myfeedforward6 as myfeedforward


#
# This program is just like tutorial3.py except a function is used to generate the html rather than a webbsite.
#
# To run 
#
#  python tutorial10.py
#
#
# This code is assuming that your python is a python 2.7 and that you have the above libraries available on your system.  
# I'm pretty sure only the requests library is not in default installation.
#
# One you get it to run, try changing the name of the mlfilename to see if you can read in different data sets.
#


#### going to read a URL ####
#mlfilename = "/Users/bylaska/bin/nn_chno_b3lyp.dat"
#mlfilename = "nn_chno_b3lyp.dat"




#### simple functions to check if string is a number ####
def evalnum(s):
   try:
      return int(s)
   except ValueError:
      return float(s)

def isevalnum(s):
   try:
      x = evalnum(s)
      return True
   except:
      return False



#### geturlresult function ####
def geturlresult(url):
   try:
      the_page = ""
      with urllib.request.urlopen(url) as response:
         the_page = response.read().rstrip()
   except:
      the_page = ""

   if isinstance(the_page,bytes): the_page = the_page.decode("utf-8")

   return the_page


covalentstr = '''
H 32 0 0 0
He 46 0 0 0
Li 133 124 0 0
Be 102 90 85 0
B 85 78 73 0
C 75 67 60 68
N 71 60 54 0
O 63 57 53 0
F 64 59 53 0
Ne 67 96 0 0
Na 155 160 0 0
Mg 139 132 127 0
Al 126 113 111 0
Si 116 107 102 0
P 111 102 94 0
S 103 94 95 0
Cl 99 95 93 0
Ar 96 107 96 0
K 196 193 0 0
Ca 171 147 133 0
Sc 148 116 114 0
Ti 136 117 108 0
V 134 112 106 0
Cr 122 111 103 0
Mn 119 105 103 0
Fe 116 109 102 0
Co 111 103 96 0
Ni 110 101 101 0
Cu 112 115 120 0
Zn 118 120 0 0
Ga 124 116 121 0
Ge 121 111 114 0
As 121 114 106 0
Se 116 107 107 0
Br 114 109 110 0
Kr 117 121 108 0
Rb 210 202 0 0
Sr 185 157 139 0
Y 163 130 124 0
Zr 154 127 121 0
Nb 147 125 116 0
Mo 138 121 113 0
Tc 128 120 110 0
Ru 125 114 103 0
Rh 125 110 106 0
Pd 120 117 112 0
Ag 128 139 137 0
Cd 136 144 0 0
In 142 136 146 0
Sn 140 130 132 0
Sb 140 133 127 0
Te 136 128 121 0
I 133 129 125 0
Xe 131 135 122 0
Cs 232 196 0 0
Ba 196 161 149 0
La 180 139 139	 0
Ce 163 137 131 0
Pr 176 138 128 0
Nd 174 137 0 0
Pm 173 135 0 0
Sm 172 134 0 0
Eu 168 134 0 0
Gd 169 135 132 0
Tb 168 135 0 0
Dy 167 133 0 0
Ho 166 133 0 0
Er 165 133 0 0
Tm 164 131 0 0
Yb 170 129 0 0
Lu 162 131 131	 0
Hf 152 128 122	 0
Ta 146 126 119	 0
W 137 120 115	 0
Re 131 119 110	 0
Os 129 116 109	 0
Ir 122 115 107	 0
Pt 123 112 110	 0
Au 124 121 123 0
Hg 133 142 0 0
Tl 144 142 150 0
Pb 144 135 137 0
Bi 151 141 135 0
Po 145 135 129 0
At 147 138 138 0
Rn 142 145 133 0
Fr 223 218 0 0
Ra 201 173 159 0
Ac 186 153 140 0
Th 175 143 136	 0
Pa 169 138 129	 0
U 170 134 118 0
Np 171 136 116 0
Pu 172 135 0  0
Am 166 135 0 0
Cm 166 136 0 0
Bk 168 139 0 0
Cf 168 140 0 0
Es 165 140 0 0
Fm 167 0 0 0
Md 173 139 0 0
No 176 0 0  0
Lr 161 141 0 0
Rf 157 140 131 0
Db 149 136 126 0
Sg 143 128 121 0
Bh 141 128 119 0
Hs 134 125 118 0
Mt 129 125 113 0
Ds 128 116 112 0
Rg 121 116 118	 0
Cn 122 137 130 0
Uut 136 0 0 0
Fl 143 0 0 0
Uup 162  0 0 0
Lv 175 0 0 0
Uus 165 0 0 0
Uuo 157  0 0 0
'''
rcovalent = {}
for ln in covalentstr.strip().split('\n'):
   ss = ln.split()
   rcovalent[ss[0]] = (0.01*eval(ss[1]),0.01*eval(ss[2]),0.01*eval(ss[3]),0.01*eval(ss[4]))


###########################################
#                                         #
#              bond_order                 #
#                                         #
###########################################
def bond_order(rc1,rc2,r12):
   dd = 0.0001
   cov = (abs(r12-(rc1[0]+rc2[0]))/(rc1[0]+rc2[0]+dd),
          abs(r12-(rc1[1]+rc2[1]))/(rc1[1]+rc2[1]+dd),
          abs(r12-(rc1[2]+rc2[2]))/(rc1[2]+rc2[2]+dd),
          abs(r12-(rc1[3]+rc2[3]))/(rc1[3]+rc2[3]+dd))
   imin = 0
   dmin = cov[0]
   if (cov[1]<dmin):
      dmin = cov[1]
      imin = 1
   if (cov[2]<dmin):
      dmin = cov[2]
      imin = 2
   if (cov[3]<dmin):
      dmin = cov[3]
      imin = 3
   b = 0
   if (cov[imin]<0.10):
      b = 1+imin
      if (imin==3):
         b = 1.5
   return b

###########################################
#                                         #
#          xyzdat_bonding_strings         #
#                                         #
###########################################

def xyzdat_bonding_strings(xyzdat):

   #### read xyz file ####
   fdict = {}
   verts  = []
   symbol = []
   rxyz   = []
   n = eval(xyzdat.strip().split("\n")[0].strip())
   for line in xyzdat.strip().split("\n")[2:]:
      if (line[1]==' '):
         key = line[0]
      else:
         key = line[0:2]
      if (key in fdict):
         fdict[key] += 1
      else:
         fdict[key] = 1
      ss = line.split()
      symbol.append(ss[0].strip())
      tple = ('atom',(0.0, ss[0].strip(), '', 0, 0, -1))
      verts.append(tple)
      rxyz.append(eval(ss[1]))
      rxyz.append(eval(ss[2]))
      rxyz.append(eval(ss[3]))

   #### generate mformula ####
   mformula = ''
   for x  in sorted(fdict.items(), key=operator.itemgetter(0)):
      mformula += x[0] + "%d" % x[1]

   #### generate adjacency matrix ####
   adjmat = []
   rij    = []
   for i in range(n):
      rij.append([0.0]*n)
      adjmat.append([0]*n)
   for i in range(n):
      for j in range(n):
         symi = symbol[i]
         symj = symbol[j]
         rci   = rcovalent[symbol[i]]
         rcj   = rcovalent[symbol[j]]
         dx = rxyz[3*i]   - rxyz[3*j]
         dy = rxyz[3*i+1] - rxyz[3*j+1]
         dz = rxyz[3*i+2] - rxyz[3*j+2]
         r = math.sqrt(dx*dx + dy*dy + dz*dz)
         rij[i][j] = r
         if i!=j:
            adjmat[i][j] = bond_order(rci,rcj,r)

   #### generate bonding ####
   covbondcount = {}
   bondcount = {}
   for i in range(n):
      for j in range(i+1,n):
         if (adjmat[i][j] > 0):
            symi = symbol[i]
            symj = symbol[j]
            if (symi<symj):
               key = symi.strip() + symj.strip()
            else:
               key = symj.strip() + symi.strip()

            if (key in bondcount):
               bondcount[key] += 1
            else:
               bondcount[key] = 1
            covkey = key + "(%.1f)" % (adjmat[i][j])
            if (covkey in covbondcount):
               covbondcount[covkey] += 1
            else:
               covbondcount[covkey] = 1
   bonding = ''
   for x  in sorted(bondcount.items(), key=operator.itemgetter(0)):
      bonding += x[0] + "%d" % x[1]
   covbonding = ''
   for x  in sorted(covbondcount.items(), key=operator.itemgetter(0)):
      covbonding += x[0] + "=%d," % x[1]
   covbonding = covbonding.strip(',')


   #### generate bonding2 ####
   bond2count = {}
   for i in range(n):
      for j in range(n):
         for k in range(j+1,n):
            if (adjmat[i][j] > 0) and (adjmat[i][k] > 0):
               symi = symbol[i]
               symj = symbol[j]
               symk = symbol[k]
               if (symj<symk):
                  key = symj.strip() + symi.strip() + symk.strip()
               else:
                  key = symk.strip() + symi.strip() + symj.strip()
               if (key in bond2count):
                  bond2count[key] += 1
               else:
                  bond2count[key] = 1
   bonding2 = ''
   for x  in sorted(bond2count.items(), key=operator.itemgetter(0)):
      bonding2 += x[0] + "%d" % x[1]

   return  (mformula + ":" + bonding + ":" + bonding2 + ":" + covbonding)






#################################
#                               #
#       read_ml_urlfile         #
#                               #
#################################
#
# This function reads an nwchem .fei file frame by frame.
# Note that urlfilename can either be a filename or a url link to a .fei file
#
def read_ml_urlfile(urlfilename):
   """Lazy function (generator) to read a fei_file 
      frame by frame."""

   if "http" in urlfilename:
      rr = requests.get(urlfilename.strip())
      mldata = rr.text.split('\n')
   else:
      with open(urlfilename,'r') as ff: 
         mldata = ff.read().strip().split('\n')

   nframes = len(mldata)
   print("nframes=",nframes)

   framecounter = 0
   while (framecounter<nframes):
      line = mldata[framecounter]
      id    = line.split()[0]
      xdata = [eval(x) for x in line.split("inputLayer:")[1].split(":inputLayer")[0].split()]
      energy= eval(line.split("outputLayer:")[1].split(":outputLayer")[0])
      framecounter += 1

      yield (id,xdata,energy)


##############################################
#                                            #
#             plot_pathenergy                #
#                                            #
##############################################
def plot_pathenergy(plot,y0,y1):
   #colors = ("blue","green","yellow","purple","orange")
   #ncc    = len(colors)
   imax      = 0 
   y0_imax   = 0.0
   y_imax    = 0.0
   ydiff_max = 0.0
   imin      = 0 
   y0_imin   = 0.0
   y_imin    = 0.0
   ydiff_min = 99e99
   #delta = 0.1
   ymin = +999999.9
   ymax = -999999.9
   for i in range(len(y1)):
      y = y1[i]
      if (y>ymax): ymax = y
      if (y<ymin): ymin = y
      if (y0[i]>ymax): ymax = y0[i]
      if (y0[i]<ymin): ymin = y0[i]

      ydiff = abs(y - y0[i])
      if (ydiff>ydiff_max):
         imax = i
         y0_imax   = y0[i]
         y_imax    = y
         ydiff_max = ydiff
      if (ydiff<ydiff_min):
         imin = i
         y0_imin   = y0[i]
         y_imin    = y
         ydiff_min = ydiff

   ### reset the plotting window and then plot data ###
   delta = 0.1*(ymax-ymin)
   plot.resetwindow(0.0,ymin-delta,1.0,ymax+delta,"Scaled Energies")

   #yb,yr = (list(t) for t in zip(*sorted(zip(y0,y1))))

   plot.plot1(y0,"red")
   plot.plot1(y1,"blue")
   #plot.dotplot1(y1,"blue")

   return ((imin,y0_imin,y_imin,ydiff_min),(imax,y0_imax,y_imax,ydiff_max))


##############################################
#                                            #
#             plot_pathenergy0               #
#                                            #
##############################################
def plot_pathenergy0(plot,y0,color):
   #colors = ("blue","green","yellow","purple","orange")
   #ncc    = len(colors)
   ymin = +999999.9
   ymax = -999999.9
   for i in range(len(y0)):
      y = y0[i]
      if (y>ymax): ymax = y
      if (y<ymin): ymin = y


   ### reset the plotting window and then plot data ###
   delta = 0.1*(ymax-ymin)
   plot.resetwindow(0.0,ymin-delta,1.0,ymax+delta,"Scaled Energies")

   #yb,yr = (list(t) for t in zip(*sorted(zip(y0,y1))))

   plot.plot1(y0,color)

   return 

##############################################
#                                            #
#             plot_energycorrelation         #
#                                            #
##############################################
def plot_energycorrelation(plot,y0,y1):
   #colors = ("blue","green","yellow","purple","orange")
   #ncc    = len(colors)
   ymin = +999999.9
   ymax = -999999.9
   for i in range(len(y1)):
      y = y1[i]
      if (y>ymax): ymax = y
      if (y<ymin): ymin = y
      if (y0[i]>ymax): ymax = y0[i]
      if (y0[i]<ymin): ymin = y0[i]

   ### reset the plotting window and then plot data ###
   delta = 0.1*(ymax-ymin)
   plot.resetwindow(ymin-delta,ymin-delta,ymax+delta,ymax+delta,"Scaled Energy Correlation")

   yb,yr = (list(t) for t in zip(*sorted(zip(y0,y1))))

   plot.dotplot(yb,yr,"black")
   #plot.dotplot1(y1,"blue")

   return 


##############################################
#                                            #
#             plot_pathdiff                  #
#                                            #
##############################################
def plot_pathdiff(plot,y0,y1):
   #ya,yb = (list(t) for t in zip(*sorted(zip(y0,y1))))

   #colors = ("blue","green","yellow","purple","orange")
   #ncc    = len(colors)
   #delta = 0.1
   ymin = +999999.9
   ymax = -999999.9
   ydiff = []
   for i in range(len(y1)):
      y =  abs(y1[i]-y0[i])
      if (y>ymax): ymax = y
      if (y<ymin): ymin = y
      ydiff.append(y)

   ### reset the plotting window and then plot data ###
   delta = 0.1*(ymax-ymin)
   plot.resetwindow(0.0,ymin,1.0,ymax+delta,"Scaled Delta Energies")
   plot.plot1(ydiff,"green")

   return (ymin,ymax)


##############################################
#                                            #
#             ytrain_generate                #
#                                            #
##############################################
def ytrain_generate(machine,weights,xdata):
   #
   ytrain = []
   for i in range(len(xdata)):
      ytrain.append(machine.evaluate(xdata[i],weights)[0])
   return ytrain


####################### main program ##############################

def main():
#
   import sys
   import getopt

   usage = \
   """
   esmiles to nn program

   Usage: mleval10 -f nn_file.dat  -n hidden layers  esmiles

   hidden_layers = "2 1"

   -h prints this message

   """

   print("\nmleval10 Arrows version")

   hidden_layers = [2,1]
   nnfilename = ""

   opts, args = getopt.getopt(sys.argv[1:], "hn:f:")
   for o, a in opts:
      if '-n' in o:
         hidden_layers = [eval(x) for x in a.strip().split()]
      if '-f' in o:
         nnfilename = a
      if o in ("-h","--help"):
         print(usage)
         exit()

   esmiles = args[0]
   print("esmiles=",esmiles)

   #weights_filename = "tutorial10"
   weights_filename = nnfilename.replace(".dat","")
   for ix in hidden_layers:
      weights_filename += "-"+str(ix) 
   weights_filename += ".weights"

   #param_filename = "tutorial10"
   param_filename = nnfilename.replace(".dat","")
   for ix in hidden_layers:
      param_filename += "-"+str(ix) 
   param_filename += ".param"


   ##### Read in params ####
   print("Param_filename:", param_filename)
   if os.path.isfile(param_filename):
      print("reading params")
      with open(param_filename,'rb') as ff:
         params = pickle.loads(ff.read())

   ninput  = params[0]
   nframes = params[1]
   emin    = params[2]
   emax    = params[3]
   emid    = params[4]
   edif    = params[5]

   print()
   print("Training Data:")
   print("ninput=",ninput)
   print("nframes=",nframes)

   print()
   print("ReScaling Energies Parameters:")
   print("emin=",emin)
   print("emax=",emax)
   print("emid=",emid)
   print("edif=",edif)


   energy_dictionary = {'C':-23725.067263640653, 'H':-320.1587671072026, 'N':-34207.75150198355, 'O':-47067.1469030725 }
   energy_type       = "gaq"

   symbol = ['C','H','N','O']
   maxatoms     = 25
   maxcharge    = 4
   maxmult      = 10

   esmiles0 = esmiles.replace(" ","%20")

   #arrows_url = 'https://arrows.emsl.pnnl.gov/api/esmiles/\"' + esmiles0.replace("/","arrowslash").strip() + '\"'
   arrows_url = 'https://arrows.emsl.pnnl.gov/api/esmiles2xyz/\"' + esmiles0.replace("/","arrowslash").strip() + '\"'

   print("nnfilename = ", nnfilename)
   print()
   print("esmiles           =",esmiles)
   print("symbols           =",symbol)
   #print("energy_type       =",energy_type)
   print("energy_dictionary =",energy_dictionary)
   print("maxatoms          =",maxatoms)
   print("maxcharge         =",maxcharge)
   print("maxmult           =",maxmult)
   print("arrows_url  =",arrows_url)
   try:
      rr = geturlresult(arrows_url)
      esmiles_dict = json.loads(rr)
   except:
      try:
         rr = geturlresult(arrows_url)
         esmiles_dict = json.loads(rr)
      except:
         print(" - API Failed")

   print("\nlen esmiles_dictionary_all=",len(esmiles_dict))
   if (len(esmiles_dict)<=0): exit()


   keys = {}
   count = 0
   n = len(symbol)

   for i in range(n):
      key = symbol[i]
      keys[key] = count
      count += 1

   for i in range(n):
      for j in range(n):
         symi = symbol[i]
         symj = symbol[j]
         if (symi<symj):
            key = symi + symj
         else:
            key = symj + symi
         if (key not in keys):
            keys[key] = count
            count += 1

   for i in range(n):
      for j in range(n):
         for k in range(n):
            symi = symbol[i]
            symj = symbol[j]
            symk = symbol[k]
            if (symj<symk):
               key = symj + symi + symk
            else:
               key = symk + symi + symj
            if (key not in keys):
               keys[key] = count
            count += 1

   print("keys=",keys)

   eoln = "\n"

   msg = eoln;
   msg += "Fetched the following entry:"
   msg += "mformula = " + esmiles_dict['mformula'] + eoln
   msg += "iupac    = " + esmiles_dict['iupac']    + eoln
   msg += "smiles   = " + esmiles_dict['smiles']   + eoln
   msg += "csmiles  = " + esmiles_dict['csmiles']  + eoln
   msg += "esmiles  = " + esmiles_dict['esmiles']  + eoln
   msg += "InChI    = " + esmiles_dict['inchi']    + eoln
   msg += "InChiKey = " + esmiles_dict['inchikey'] + eoln
   #msg += "cid      = " + esmiles_dict['cid']      + eoln
   #msg += "cas      = " + str(esmiles_dict['cas'])      + eoln
   #if (esmiles_dict['kegg'] is not None): msg += "kegg     = " + esmiles_dict['kegg']     + eoln
   msg += "bonding_string  = " + esmiles_dict['bonding_string']  + eoln
   msg += "covalent_string = " + esmiles_dict['covalent_string'] + eoln
   msg += "charge          = " + str(esmiles_dict['charge'])     + eoln
   msg += "mult            = " + str(esmiles_dict['mult'])       + eoln
   print(msg)

   charge   = esmiles_dict['charge']
   mult     = esmiles_dict['mult']
   xyz_blob = esmiles_dict['xyz_blob']

   msg2 = xyzdat_bonding_strings(xyz_blob)
   print("msg2=",msg2)
   bonding1 = re.split('(\d+)', msg2.strip().split(":")[0])
   bonding2 = re.split('(\d+)', msg2.strip().split(":")[1])
   bonding3 = re.split('(\d+)', msg2.strip().split(":")[2])

   print("Bonding1=",bonding1)
   print("Bonding2=",bonding2)
   print("Bonding3=",bonding3)

   nkeys = len(keys)
   input_hash = [0]* (nkeys+ (2*maxcharge+1) + maxmult)
   input_hash[nkeys + maxcharge + charge] = 1
   input_hash[nkeys + (2*maxcharge+1) + (mult-1)] = 1
   x1 = bonding1[0:-1:2] + bonding2[0:-1:2] + bonding3[0:-1:2]
   y1 = bonding1[1::2]   + bonding2[1::2]   + bonding3[1::2]

   keystop = False
   for x in x1:
      if (x not in keys): keystop = True
   if keystop:
       print("Exiting because key not in symbols\n")
       exit()

   for i in range(len(x1)):
      input_hash[keys[x1[i]]] = eval(y1[i])

   eref = 0.0
   for sym in symbol:
      nsym = input_hash[keys[sym]]
      eref += nsym*energy_dictionary[sym]

   natoms       = eval(xyz_blob.strip().split("\n")[0])
   inputlayer   = [0]*len(input_hash)
   for i in range(len(inputlayer)): inputlayer[i] = input_hash[i]
   for i in range(nkeys):           inputlayer[i] = input_hash[i]/natoms

   #outputlayer         = (output_hash - eref)/natoms/627.509469

   if (natoms>maxatoms):
       print("Exiting because natoms>maxatoms\n")
       exit()

   print("input_hash =",input_hash)
   #print("output_hash=",output_hash)

   print()
   print()
   print("NN Input:")
   #print("Id                =",Id)
   #print("esmiles           =",outesmiles)
   #print("theory            =",theory)
   #print("xc                =",xc)
   #print("basis             =",basis)
   #print("charge            =",charge)
   #print("solvation_type    =",solvation_type)
   print("mult              =",mult)
   print("natoms            =",natoms)
   for sym in symbol:
      print("number of " + sym + "-atoms = ",input_hash[keys[sym]],energy_dictionary[sym])
   print("eref (kcal/mol)   =",eref)
   print("ML E to kcal/mol  =",natoms*627.509469)

   #print("eref              =",eref,(output_hash-eref),(output_hash-eref)/natoms/627.509469)

   print()
   print("nkeys             =",len(keys))
   print("keys              =",keys)
   print()
   print("sizeqm            =",2*maxcharge+1+maxmult)
   print("len(inputlayer)   =",len(inputlayer))
   print()
   print("InputLayer        =",inputlayer)
   #print("OutputLayer       =",outputlayer)
   print()

   nnstr = " inputLayer: "
   for x in inputlayer:
      nnstr += str(x) + " "
   nnstr += ":inputLayer "
   print(nnstr)

   #beta = 2.0
   #sigmoid   = lambda x: 1.0/(1.0+math.exp(-x))
   #sigmoidp  = lambda x: math.exp(-x)/(1.0+math.exp(-x))**2
   #sigmoidpp = lambda x: math.exp(-x)*(math.exp(-x)-1.0)/(1.0+math.exp(-x))**3
   #ap = 1.0
   #xp = 4.5
   #bp = 3.0
   #penalty  = lambda x: ap*(0.5*(math.tanh(bp*(x-xp)) - math.tanh(bp*(x+xp))) + 1.0)
   #penaltyp = lambda x: ap*0.5*bp*( (1/math.cosh(bp*(x-xp)))**2 - (1.0/math.cosh(bp*(x+xp)))**2)
   #
   #sigmoid   = lambda x: 0.5*(math.tanh(beta*x)+1.0)
   #sigmoidp  = lambda x: 0.5*beta*(1.0/math.cosh(beta*x))**2
   #sigmoidpp = lambda x: 0.5*(-2.0)*beta*beta*math.tanh(beta*x)*(1.0/math.sech(beta*x))**2
 
   xmoid1    = lambda x: x
   xmoidp1   = lambda x: 1.0
   xmoidpp1  = lambda x: 0.0

   relu      = lambda x: max(0.0,x)
   relup     = lambda x: 0.0 if (x<=0.0) else 1.0
   relupp    = lambda x: 0.0

   #bias = [[0.01],[0.01],[0.001],[0.0001],[0.00001],[0.0000001]]
   #bias = [[0.01],[0.01],[0.001]]
   #bias = []

   ### check/create network topology ###
   print()
   print("creating network topology - ninput x hidden_layers x 1")
   print("                          - hidden_layers = ( ninput *",hidden_layers,")")
   print()

   ### define the network topology ###
   ntwrk0 = [ninput]
   ntwrkf = [xmoid1]
   ntwrkp = [xmoidp1]
   ntwrkpp= [xmoidpp1]
   for h in range(len(hidden_layers)):
      ntwrk0.append(hidden_layers[h]*ninput)
      ntwrkf.append(relu)
      ntwrkp.append(relup)
      ntwrkpp.append(relupp)
   ntwrk0.append(1)
   ntwrkf.append(xmoid1)
   ntwrkp.append(xmoidp1)
   ntwrkpp.append(xmoidpp1)

   #network = [[ninput,2*ninput,ninput,1],[xmoid1,relu,relu,xmoid1],[xmoidp1,relup,relup,xmoidp1],[xmoidpp1,relupp,relupp,xmoidpp1]]
   #machine = myfeedforward.MyFeedForward(ninput,network[0],network[1],network[2],network[3])

   machine = myfeedforward.MyFeedForward(ntwrk0,ntwrkf,ntwrkp,ntwrkpp)

   print("Weights_filename:", weights_filename)
   if os.path.isfile(weights_filename):
      print("reading weights")
      with open(weights_filename,'rb') as ff:
         weights = pickle.loads(ff.read())
      

   nw = len(weights)
   print("nw=",nw)


   ### generate current ytrain ###
   ytrain=machine.evaluate(inputlayer,weights)

   print("ytrain=",ytrain)

   scaled_energy = ytrain[0]
   energy = 2*edif*scaled_energy + emin
   final_energy = natoms*627.509469*energy + eref


   print("energy=",energy)
   print("final energy=",final_energy,final_energy/23.06/27.2116)

   return


if __name__ == "__main__":
   main()


