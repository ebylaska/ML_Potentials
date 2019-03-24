from __future__ import print_function

### Need to make sure these libraries can be loaded! ###
import os,sys,subprocess,urllib2,requests,getopt,pickle
import math,random
import xyplotter,webbrowser
import myfeedforward5 as myfeedforward


#
# This program is just like tutorial3.py except a function is used to generate the html rather than a webbsite.
#
# To run 
#
#  python tutorial5.py
#
#
# This code is assuming that your python is a python 2.7 and that you have the above libraries available on your system.  
# I'm pretty sure only the requests library is not in default installation.
#
# One you get it to run, try changing the name of the feifilename to see if you can read in different data sets.
#


#feifilename = "h2-aimd.fei"
#feifilename = "tequil-2018-3-2-12.fei"
#feifilename = "tequil-2018-3-9-8.fei"

#### going to read a URL ####
#feifilename = "https://arrows.emsl.pnnl.gov/api/eric_view/raw=we31869:/media/Seagate2/Projects/ForJim/Position-Specific-Isotopes/1M/Pyruvate/AIMD/tequil-2018-3-2-12.fei"
#feifilename = "https://arrows.emsl.pnnl.gov/api/eric_view/raw=we31869:/media/Seagate2/Projects/BES/Mackinawite/Cascade-hopper/udimer-fes.fei"
feifilename = "https://arrows.emsl.pnnl.gov/api/eric_view/raw=we31869:/media/Seagate2/Projects/ForRaymond/perm/co2-2019-30-19-12.fei"




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
#       read_fei_urlfile        #
#                               #
#################################
#
# This function reads an nwchem .fei file frame by frame.
# Note that urlfilename can either be a filename or a url link to a .fei file
#
def read_fei_urlfile(urlfilename):
   """Lazy function (generator) to read a fei_file 
      frame by frame."""
   #tobohr = 1.0/0.529177
   tobohr = 1.0  #fei file is in atomic units

   if "http" in urlfilename:
      rr = requests.get(urlfilename.strip())
      xyzdata = rr.text.split('\n')
   else:
      with open(filename,'r') as ff: 
         xyzdata = ff.read().split('\n')

   nion = int(xyzdata[0])
   nlines = nion+5
   nframes = len(xyzdata)/nlines
   print("NION=",nion," nframes=",nframes)

   framecounter = 0
   while (framecounter<nframes):
      lines   = xyzdata[framecounter*nlines:(framecounter+1)*nlines]

      energy  = float(lines[1])
      symbols = [ln.split()[1] for ln in lines[5:]]
      rions   = [tobohr*float(ln.split()[i]) for ln in lines[5:] for i in range(3,6)]
      fions   = [tobohr*float(ln.split()[i]) for ln in lines[5:] for i in range(6,9)]
      framecounter += 1

      yield (symbols,rions,fions,energy)


##############################################
#                                            #
#             plot_pathenergy                #
#                                            #
##############################################
def plot_pathenergy(plot,y0,machine,weights,ds1,ds2,ds3):
   #colors = ("blue","green","yellow","purple","orange")
   #ncc    = len(colors)
   delta = 0.1
   ymin = +999999.9
   ymax = -999999.9
   y1 = []
   for i in range(len(ds1)):
      y = machine.evaluate([ds1[i],ds2[i],ds3[i]],weights)[0]
      if (y>ymax): ymax = y
      if (y<ymin): ymin = y
      if (y0[i]>ymax): ymax = y0[i]
      if (y0[i]<ymin): ymin = y0[i]
      y1.append(y)

   ### reset the plotting window and then plot data ###
   plot.resetwindow(0.0,ymin-delta,1.0,ymax+delta,"Scaled Energies")
   plot.plot1(y0,"blue")
   plot.plot1(y1,"red")
   #plot.dotplot1(y1,"blue")




##############################################
#                                            #
#             xydata_plotdatajs              #
#                                            #
##############################################
# This function generates html code to plot
# data.
def xydata_plotdatajs(edat,label):

   title = label
   xylabels = ''
   if ('#Title' in edat):
      title = edat.split("#Title")[1].split('\n')[0]
   else:
      title = label
   if ('#Labels' in edat): xylabels = edat.split("#Labels")[1].split('\n')[0]

   elist = edat.strip().split("\n")
   while ("#" in elist[0]):
      elist = elist[1:]

   ny = len(elist[0].split())-1
   hasheader = not isevalnum(elist[0].split()[0])

   if (hasheader):
      xlabel = elist[0].split()[0]
      xdat = "['%s'," % (elist[0].split()[0])
      ydat = []
      for i in range(ny):
         ydat.append("['%s'," % (elist[0].split()[i+1]))
      elist = elist[1:]
   else:
      xlabel = 'x'
      xdat = "['x',"
      ydat = []
      for i in range(ny):
         ydat.append("['y%d'," % i)
      if (xylabels!=''):
         xylist = xylabels.split()
         xlabel = xylist[0]
         xdat = "['%s'," % (xylist[0])
         for i in range(1,len(xylist)):
            ydat[i-1] = "['%s'," % (xylist[i])


   for ee in elist:
      ss = ee.split()
      if ("#" not in ss[0]):
         xdat += ss[0] + ", "
         for i in range(ny):
            ydat[i] += ss[i+1] + ", "
   xdat  = xdat.rstrip(',') + "]"
   for i in range(ny):
      ydat[i] = ydat[i].rstrip(',') + "]"

   msg4 = "<html>\n"
   msg4 += '''
    <link href="https://cdnjs.cloudflare.com/ajax/libs/c3/0.6.9/c3.min.css" rel="stylesheet">
    <script type="text/javascript" src="https://d3js.org/d3.v5.min.js" charset="utf-8"></script>
    <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/c3/0.6.9/c3.min.js"></script>
    <br> <center><b> %s </b></center>
    <div id="chart"></div>
    <script type="text/javascript">
    var chart = c3.generate({
       bindto: '#chart',
       size: { height: 480},
       data: {
         type: 'spline',
         x: '%s',
         columns: [
           %s  ''' % (title,xlabel,xdat)
   for i in range(ny): msg4 += ", %s" % (ydat[i])
   msg4 += '''
         ]
       },
       axis: {
           x: {
               label: '%s',
               tick: {count: 10, format: d3.format(".2f"), culling: false}
           },
           y: {
               label: 'y data'
           }
       }
    });
    </script>
    ''' % (xlabel)

   msg4 += "</html>\n"

   return msg4




####################### main program ##############################

### started from scratch - build up path from l=0 ###
plot  = xyplotter.xyplotter(0.0,0.0,1.0,1.0, "Scaled Energies Plot",2)


##### Read in fei file lazilly and do a simple print out of the data ####
nframes = 0

emin = +99999.99
emax = -99999.99
energies = []

d1min = +999999.9
d1max = -999999.9
dist1 = []

d2min = +999999.9
d2max = -999999.9
dist2 = []

d3min = +999999.9
d3max = -999999.9
dist3 = []
for (symbols,rions,fions,energy) in read_fei_urlfile(feifilename):
   energies.append(energy)
   if (energy<emin): emin = energy
   if (energy>emax): emax = energy

   x12 = rions[0] - rions[3]
   y12 = rions[1] - rions[4]
   z12 = rions[2] - rions[5]
   r12 = math.sqrt(x12*x12 + y12*y12 + z12*z12)
   dist1.append(r12)
   if (r12<d1min): d1min = r12
   if (r12>d1max): d1max = r12

   x23 = rions[3] - rions[6]
   y23 = rions[4] - rions[7]
   z23 = rions[5] - rions[8]
   r23 = math.sqrt(x23*x23 + y23*y23 + z23*z23)
   dist2.append(r23)
   if (r23<d2min): d2min = r23
   if (r23>d2max): d2max = r23

   x13 = rions[0] - rions[6]
   y13 = rions[1] - rions[7]
   z13 = rions[2] - rions[8]
   r13 = math.sqrt(x13*x13 + y13*y13 + z13*z13)
   dist3.append(r13)
   if (r13<d3min): d3min = r13
   if (r13>d3max): d3max = r13

   nframes += 1

#es = (e-emin)/(emax-emin)
#x = ((xmax+xmin) + xs*(xmax-xmin))*0.5
#xs = (2*x - (xmax+xmin))/(xmax-xmin)
#
#define scaled coordinates ###
for i in range(nframes):
   energies[i] = (energies[i]-emin)/(emax-emin)
   dist1[i] = (2*dist1[i] - (d1max+d1min))/(d1max-d1min)
   dist2[i] = (2*dist2[i] - (d2max+d2min))/(d2max-d2min)
   dist3[i] = (2*dist3[i] - (d3max+d3min))/(d3max-d3min)

nframes0 = nframes-500


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

sigmoid      = lambda x: 0.5*(math.tanh(beta*x)+1.0)
sigmoidp     = lambda x: 0.5*beta*(1.0/math.cosh(beta*x))**2
sigmoidpp    = lambda x: 0.5*(-2.0)*beta*beta*math.tanh(beta*x)*(1.0/math.sech(beta*x))**2
xmoid1   = lambda x: x
xmoidp1  = lambda x: 1.0
xmoidpp1 = lambda x: 0.0

#bias = [[0.01],[0.01],[0.001],[0.0001],[0.00001],[0.0000001]]
bias = [[0.01],[0.01],[0.001]]
bias = []
machine = myfeedforward.MyFeedForward([3,60,120,1],[xmoid1,sigmoid,sigmoid,xmoid1],[xmoidp1,sigmoidp,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,sigmoidpp,xmoidpp1],bias)

if os.path.isfile("tutorial5.weights"):
   with open("tutorial5.weights",'r') as ff:
      weights = pickle.loads(ff.read())
else:
   weights = machine.initial_w()
   for i in range(len(weights)):
      weights[i] *= 1.0

nw = len(weights)

print("weights=",weights)
m = [0.0]*len(weights)
v = [0.0]*len(weights)
beta1t = 1.0
beta2t = 1.0
print("nframes=",nframes)



###plot the relative energies using Turtle graphics###
plot_pathenergy(plot,energies[nframes0:],machine,weights,dist1[nframes0:],dist2[nframes0:],dist3[nframes0:])

enter0 = raw_input(" -- start simulation --- ")

error = 999999.9
ii = 0
for i in range(1000000):

   es = energies[ii]
   xs = dist1[ii]
   ys = dist2[ii]
   us = dist3[ii]
   ii = ((ii+1)%nframes0)
   
   gg = machine.w_energy_gradient([xs,ys,us],[es],weights)
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

   if ((i%1000)==0):
      es1 = machine.evaluate([xs,ys,us],weights)[0]
      print("%10d %10.3f %10.3f %10.3f %10.3f %10.3f  ||  %12.6f" % (i,xs,ys,us,es,es1,error))
      plot_pathenergy(plot,energies[nframes0:],machine,weights,dist1[nframes0:],dist2[nframes0:],dist3[nframes0:])

      with open("tutorial5.weights",'w') as ff:
         ff.write(pickle.dumps(weights))

###plot the relative energies using html ###

###wait for return so that plot can be seen###
x = raw_input("--Press return to finish--")
print("yeh - writing weeights")

with open("tutorial5.weights",'w') as ff:
   ff.write(pickle.dumps(weights))

