from __future__ import print_function

### Need to make sure these libraries can be loaded! ###
import os,sys,subprocess,urllib2,requests,getopt
import xyplotter


#
# This program extends upon tutorial1.py to demonstrate how to use a simple python plotter written with turtle graphics.
#
# To run 
#
#  python tutorial2.py
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
feifilename = "https://arrows.emsl.pnnl.gov/api/eric_view/raw=we31869:/media/Seagate2/Projects/BES/Mackinawite/Cascade-hopper/udimer-fes.fei"



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
def plot_pathenergy(plot,energies):
   #colors = ("blue","green","yellow","purple","orange")
   #ncc    = len(colors)
   delta = 0.0
   y = []
   y0 = energies[0]
   ymax = -999999.9e+99
   ymin = +999999.9e+99
   for e in energies:
      y1 = (e - y0)*27.2116*23.06
      if (y1>ymax): ymax = y1
      if (y1<ymin): ymin = y1
      y.append(y1)

   ### reset the plotting window and then plot data ###
   plot.resetwindow(0.0,ymin-delta,1.0,ymax+delta,"Relative Energies")
   plot.plot1(y,"blue")



####################### main program ##############################

### started from scratch - build up path from l=0 ###
plot  = xyplotter.xyplotter(0.0,0.0,1.0,1.0, "Simple Plot for Raymond",2)


##### Read in fei file lazilly and do a simple print out of the data ####
frame = 0
energies = []
for (symbols,rions,fions,energy) in read_fei_urlfile(feifilename):
   energies.append(energy)
   #nion = len(symbols)
   #print("frame=",frame)
   #print("energy=",energy)
   #print("nion=",nion)
   #print("symbols=",symbols)
   #print("rions=",rions)
   #print("fions=",rions)
   frame += 1


###plot the relative energies ###
plot_pathenergy(plot,energies)


###wait for return so that plot can be seen###
x = raw_input("--Press return to finish--")
print("yeh")

