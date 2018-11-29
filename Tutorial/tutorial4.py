from __future__ import print_function

### Need to make sure these libraries can be loaded! ###
import os,sys,subprocess,urllib2,requests,getopt
import xyplotter,webbrowser


#
# This program is just like tutorial3.py except a function is used to generate the html rather than a webbsite.
#
# To run 
#
#  python tutorial4.py
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


###plot the relative energies using Turtle graphics###
plot_pathenergy(plot,energies)


###plot the relative energies using html ###
data = "#Title A simpler html plot for Raymond\n"
data += "#Labels frame_number Energy\n"
for i in range(len(energies)):
   data += "%f %f\n" % (i*1.0,energies[i])

html = xydata_plotdatajs(data,"XY plot")

path = os.path.abspath('tutorial4.html')
url = 'file://' + path
with open(path, 'w') as ff:
   ff.write(html)
webbrowser.open(url)


###wait for return so that plot can be seen###
x = raw_input("--Press return to finish--")
print("yeh")

