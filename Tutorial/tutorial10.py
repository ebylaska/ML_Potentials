
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

   #delta = 0.1
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
   delta = 0.1*(ymax-ymin)
   plot.resetwindow(0.0,ymin-delta,1.0,ymax+delta,"Scaled Energies")
   plot.plot1(y0,"blue")
   plot.plot1(y1,"red")
   #plot.dotplot1(y1,"blue")

   return ((imin,y0_imin,y_imin,ydiff_min),(imax,y0_imax,y_imax,ydiff_max))


##############################################
#                                            #
#             plot_pathdiff                  #
#                                            #
##############################################
def plot_pathdiff(plot,y0,machine,weights,xdata):
   #colors = ("blue","green","yellow","purple","orange")
   #ncc    = len(colors)

   #delta = 0.1
   ymin = +999999.9
   ymax = -999999.9
   y1 = []
   for i in range(len(xdata)):
      y = abs(machine.evaluate(xdata[i],weights)[0]-y0[i])
      if (y>ymax): ymax = y
      if (y<ymin): ymin = y
      y1.append(y)

   ### reset the plotting window and then plot data ###
   delta = 0.1*(ymax-ymin)
   plot.resetwindow(0.0,ymin,1.0,ymax+delta,"Scaled Delta Energies")
   plot.plot1(y1,"green")

   return (ymin,ymax)



####################### main program ##############################

def main():
#
   import sys
   import getopt

   usage = \
   """
   esmiles to nn program

   Usage: tutorial10 -n "hidden layers"

   hidden_layers = "2 1"

   -h prints this message

   """

   print("\ntutorial10 Arrows version")

   hidden_layers = [2,1]

   opts, args = getopt.getopt(sys.argv[1:], "hn:")
   for o, a in opts:
      if '-n' in o:
         hidden_layers = [eval(x) for x in a.strip().split()]
      if o in ("-h","--help"):
         print(usage)
         exit()

   weights_filename = "tutorial10"
   for ix in hidden_layers:
      weights_filename += "-"+str(ix) 
   weights_filename += ".weights"


   ##### Read in fei file lazilly and do a simple print out of the data ####
   nframes = 0
   emin = +99999.99
   emax = -99999.99
   energies = []
   xinputs  = []
   ids      = []
   for (id,xdata,energy) in read_ml_urlfile(mlfilename):
      #print("nframes=",nframes)
      xinputs.append(xdata)
      energies.append(energy)
      ids.append(id)
      if (energy<emin): emin = energy
      if (energy>emax): emax = energy

      nframes += 1

   emid = (emax+emin)/2.0
   edif = (emax-emin)/2.0

   ninput   = len(xinputs[0])
   nframes0 = 0

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

   scaled_energies = [0.0]*nframes
   for i in range(nframes):
      scaled_energies[i] = (energies[i]-emin)/(2*edif)
    


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
   else:
      print("creating initial weights")
      weights = machine.initial_w()
      for i in range(len(weights)):
         weights[i] *= 1.0

   nw = len(weights)
   print("nw=",nw)


   nepoch = 100
   nbatch = 25
   print()
   print("Minimization Parameters:")
   print("nepoch =",nepoch)
   print("nbatch =",nbatch)


   ###plot the relative energies using Turtle graphics###
   plot  = xyplotter.xyplotter(0.0,0.0,1.0,1.0, "Scaled Energies Plot",2)
   iydiff = plot_pathenergy(plot,scaled_energies[nframes0:],machine,weights,xinputs[nframes0:])
   id_min = ids[iydiff[0][0]]
   id_max = ids[iydiff[1][0]]
   print("iydiff=",iydiff," idmin=",id_min, " idmax=",id_max)
   plot2 = xyplotter.xyplotter(0.0,0.0,1.0,1.0, "Scaled Diff Energies Plot",2)
   yminmax = plot_pathdiff(plot2,scaled_energies[nframes0:],machine,weights,xinputs[nframes0:])

   enter0 = input(" -- start simulation --- ")

   alpha = 0.0001
   beta1 = 0.9
   beta2 = 0.999
   eps   = 1e-8
   beta1t = 1.0
   beta2t = 1.0
   m    = [0.0]*nw
   v    = [0.0]*nw
   gsum = [0.0]*nw

   for epoch in range(nepoch):
      batch   = 0
      error1 = 0.0
      for j in range(nw): gsum[j] = 0.0
   
      for i in range(nframes):
         batch += 1
         gg = machine.w_energy_gradient(xinputs[i],[scaled_energies[i]],weights)
         error = gg[0]
         g1    = gg[1]

         error1 += error
         for j in range(nw): 
            gsum[j] += g1[j]

         if ((batch>=nbatch) or (i>=(nframes-1))):
            xx = 1.0/float(batch)
            for j in range(nw): 
               gsum[j] *= xx

            beta1t *= beta1
            beta2t *= beta2
            if (beta1t<1.0e-20): beta1t = 0.0
            if (beta2t<1.0e-20): beta2t = 0.0
            alphat = alpha*math.sqrt(1.0-beta2t)/(1.0-beta1t)

            for j in range(nw):
               m[j] = beta1*m[j] + (1.0-beta1)*gsum[j]
               v[j] = beta2*v[j] + (1.0-beta2)*(gsum[j]*gsum[j])
               weights[j] -= alphat*m[j]/(math.sqrt(v[j]) + eps)

            batch = 0
            for j in range(nw): gsum[j] = 0.0

      ### plot the relative scaled_energies using turtle graphicsl ###
      iydiff = plot_pathenergy(plot,scaled_energies[nframes0:],machine,weights,xinputs[nframes0:])
      id_min = ids[iydiff[0][0]]
      id_max = ids[iydiff[1][0]]
      yminmax = plot_pathdiff(plot2,scaled_energies[nframes0:],machine,weights,xinputs[nframes0:])

      ### print out update ###
      msg  = "epoch=%5d nframes=%5d nbatch=%3d " % (epoch+1,nframes,nbatch)
      msg += "| min=%5d %6.3f %6.3f %8.3e %5s " % (iydiff[0][0],iydiff[0][1],iydiff[0][2],iydiff[0][3],id_min)
      msg += "| max=%5d %6.3f %6.3f %8.3e %5s " % (iydiff[1][0],iydiff[1][1],iydiff[1][2],iydiff[1][3],id_max)
      msg += "|| error = %12.8e" % (math.sqrt(error1))
      print(msg)

      ### writeout weights ###
      with open(weights_filename,'wb') as ff:
         ff.write(pickle.dumps(weights))

   ### end for epoch ###


   ### wait for return so that plot can be seen ###
   x = input("--Press return to finish--")

   ### writeout weights ###
   with open(weights_filename,'wb') as ff:
      ff.write(pickle.dumps(weights))

   return


if __name__ == "__main__":
   main()


