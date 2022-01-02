from __future__ import print_function

### Need to make sure these libraries can be loaded! ###
import os,sys,subprocess,urllib,requests,getopt,pickle
import math,random
import xyplotter,webbrowser
import myfeedforward6 as myfeedforward


#
# This program is just like tutorial3.py except a function is used to generate the html rather than a webbsite.
#
# To run 
#
#  python tutorial9.py
#
#
# This code is assuming that your python is a python 2.7 and that you have the above libraries available on your system.  
# I'm pretty sure only the requests library is not in default installation.
#
# One you get it to run, try changing the name of the mlfilename to see if you can read in different data sets.
#


#### going to read a URL ####
mlfilename = "/Users/bylaska/bin/nn_chno_b3lyp.dat"
mlfilename = "nn_chno_b3lyp.dat"




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
def plot_pathenergy(plot,y0,machine,weights,xdata):
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

   delta = 0.1
   ymin = +999999.9
   ymax = -999999.9
   y1 = []
   for i in range(len(xdata)):
      y = machine.evaluate(xdata[i],weights)[0]
      if (y>ymax): ymax = y
      if (y<ymin): ymin = y
      if (y0[i]>ymax): ymax = y0[i]
      if (y0[i]<ymin): ymin = y0[i]
      y1.append(y)

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
   plot.resetwindow(0.0,ymin-delta,1.0,ymax+delta,"Scaled Energies")
   plot.plot1(y0,"blue")
   plot.plot1(y1,"red")
   #plot.dotplot1(y1,"blue")

   return ((imin,y0_imin,y_imin,ydiff_min),(imax,y0_imax,y_imax,ydiff_max))





####################### main program ##############################

### started from scratch - build up path from l=0 ###


##### Read in fei file lazilly and do a simple print out of the data ####
nframes = 0

emin = +99999.99
emax = -99999.99
energies = []
xinputs  = []
ids      = []
for (id,xdata,energy) in read_ml_urlfile(mlfilename):
   print("nframes=",nframes)
   xinputs.append(xdata)
   energies.append(energy)
   ids.append(id)
   if (energy<emin): emin = energy
   if (energy>emax): emax = energy

   nframes += 1

nframes0 = 0


alpha0 = 0.01
alpha = 0.0001
beta1 = 0.9
beta2 = 0.999
eps   = 1e-8

beta = 2.0
#sigmoid   = lambda x: 1.0/(1.0+math.exp(-x))
#sigmoidp  = lambda x: math.exp(-x)/(1.0+math.exp(-x))**2
#sigmoidpp = lambda x: math.exp(-x)*(math.exp(-x)-1.0)/(1.0+math.exp(-x))**3
ap = 1.0
xp = 4.5
bp = 3.0
penalty  = lambda x: ap*(0.5*(math.tanh(bp*(x-xp)) - math.tanh(bp*(x+xp))) + 1.0)
penaltyp = lambda x: ap*0.5*bp*( (1/math.cosh(bp*(x-xp)))**2 - (1.0/math.cosh(bp*(x+xp)))**2)


sigmoid   = lambda x: 0.5*(math.tanh(beta*x)+1.0)
sigmoidp  = lambda x: 0.5*beta*(1.0/math.cosh(beta*x))**2
sigmoidpp = lambda x: 0.5*(-2.0)*beta*beta*math.tanh(beta*x)*(1.0/math.sech(beta*x))**2

xmoid1    = lambda x: x
xmoidp1   = lambda x: 1.0
xmoidpp1  = lambda x: 0.0

relu      = lambda x: max(0.0,x)
relup     = lambda x: 0.0 if (x<=0.0) else 1.0
relupp    = lambda x: 0.0

#bias = [[0.01],[0.01],[0.001],[0.0001],[0.00001],[0.0000001]]
bias = [[0.01],[0.01],[0.001]]
bias = []

ninput = len(xinputs[0])
print("ninput=",ninput)
machine = myfeedforward.MyFeedForward([ninput,2*ninput,2*ninput,1],[xmoid1,relu,relu,xmoid1],[xmoidp1,relup,relup,xmoidp1],[xmoidpp1,relupp,relupp,xmoidpp1])

if os.path.isfile("tutorial9.weights"):
   with open("tutorial9.weights",'rb') as ff:
      weights = pickle.loads(ff.read())
else:
   weights = machine.initial_w()
   for i in range(len(weights)):
      weights[i] *= 1.0

nw = len(weights)

print("weights=",weights)
print("nw=",nw)
m = [0.0]*len(weights)
v = [0.0]*len(weights)
beta1t = 1.0
beta2t = 1.0
print("nframes=",nframes)



y = machine.evaluate(xinputs[0],weights)[0]
print("y=",y)

gg = machine.w_energy_gradient(xinputs[0],[energies[0]],weights)
print("error=",gg[0])

###plot the relative energies using Turtle graphics###
plot  = xyplotter.xyplotter(0.0,0.0,1.0,1.0, "Scaled Energies Plot",2)
iydiff = plot_pathenergy(plot,energies[nframes0:],machine,weights,xinputs[nframes0:])
id_min = ids[iydiff[0][0]]
id_max = ids[iydiff[1][0]]
print("iydiff=",iydiff," idmin=",id_min, " idmax=",id_max)

enter0 = input(" -- start simulation --- ")

error = 999999.9
ii = 0
for i in range(1000000):

   es = energies[ii]
   ii = ((ii+1)%nframes)
   
   gg = machine.w_energy_gradient(xinputs[ii],[energies[ii]],weights)
   error0 = error
   error = gg[0]
   g1    = gg[1]

   #for j in range(nw):
   #   error += penalty(weights[j])
   #   g1[j] += penaltyp(weights[j])


   beta1t *= beta1
   beta2t *= beta2
   alphat = alpha*math.sqrt(1.0-beta2t)/(1.0-beta1t)
   for j in range(nw):
      m[j] = beta1*m[j] + (1.0-beta1)*g1[j]
      v[j] = beta2*v[j] + (1.0-beta2)*g1[j]*g1[j]
      weights[j] -= alphat*m[j]/(math.sqrt(v[j]) + eps)

   if ((i%5000)==0):
      es1 = machine.evaluate(xinputs[ii],weights)[0]
      iydiff = plot_pathenergy(plot,energies[nframes0:],machine,weights,xinputs[nframes0:])
      id_min = ids[iydiff[0][0]]
      id_max = ids[iydiff[1][0]]
      print("%10d cur=%5d %6.3f %6.3f %8.3e %5s  min=%5d %6.3f %6.3f %8.3e %5s  max=%5d %6.3f %6.3f %8.3e %5s ||  error = %12.8e" % (i,ii,es,es1,abs(es-es1),ids[ii],iydiff[0][0],iydiff[0][1],iydiff[0][2],iydiff[0][3],id_min,iydiff[1][0],iydiff[1][1],iydiff[1][2],iydiff[1][3],id_max,math.sqrt(error)))

      with open("tutorial9.weights",'wb') as ff:
         ff.write(pickle.dumps(weights))

###plot the relative energies using html ###

###wait for return so that plot can be seen###
x = input("--Press return to finish--")
print("yeh - writing weeights")

with open("tutorial9.weights",'wb') as ff:
   ff.write(pickle.dumps(weights))

