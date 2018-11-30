from __future__ import print_function

import os,sys,subprocess,urllib2,requests,getopt,webbrowser,math,xyplotter
import fastatomfm3 as atomfm

#
# This program demonstrates how to genreate feature functions
#
# To run 
#
#  python fei2ff.py
#
#
# This code is assuming that your python is a python 2.7 and that you have the above libraries available on your system.  
# I'm pretty sure only the requests library is not in default installation.
#
# One you get it to run, try changing the name of the feifilename to see if you can read in different data sets.
#


#xyzfilename = "small.xyz"
#xyzfilename = "../eatomstructure.xyz"
#feifilename = "h2-aimd.fei"
#feifilename = "tequil-2018-3-2-12.fei"
feifilename = "tequil-2018-3-9-8.fei"


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

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
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
#             plot_ffpaths                   #
#                                            #
##############################################
def plot_ffpaths(plot,data):
   colors = ("blue","green","yellow","purple","orange")
   ncc    = len(colors)
   delta = 0.0
   yplot = None
   ymax = -999999.9e+99
   ymin = +999999.9e+99
   for line in data.split('\n'):
      ss = line.split()
      if (len(ss)>1):
         if (yplot is None):
            yplot = []
            for s in ss[1:]:
               y1 = eval(s) 
               if (y1>ymax): ymax = y1
               if (y1<ymin): ymin = y1
               yplot.append([y1])
         else:
            ii = 0
            for s in ss[1:]:
               y1 = eval(s) 
               if (y1>ymax): ymax = y1
               if (y1<ymin): ymin = y1
               yplot[ii] += [y1]
               ii += 1
   plot.resetwindow(0.0,ymin-delta,1.0,ymax+delta,"Plot of feature functions")
   pc=0; pcc=0;
   for yy in yplot: 
      plot.plot1(yy,colors[pc])
      pcc += 1
      if ((pcc%ncc)==0): pc = (pc+1)%ncc




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





#### read in parameters ####
parameters = []
with open("parameters",'r') as ff: 
   data = ff.read()
   for line in data.split('\n'):
      ss = line.split()
      if (len(ss)>1):
         p = []
         for s in ss:
            if s.isdigit():  
               p += [int(s)]
            elif isFloat(s): 
               p += [float(s)]
         ok = p[0] in range(1,6)
         if (p[0]==1): ok = ok and (len(p)==2)
         if (p[0]==2): ok = ok and (len(p)==4)
         if (p[0]==3): ok = ok and (len(p)==3)
         if (p[0]==4): ok = ok and (len(p)==5)
         if (p[0]==5): ok = ok and (len(p)==5)
         if ok: parameters.append(p)

print("#parameters=",parameters,len(parameters))

#### define Atomic Feature Mapping ####
afm = atomfm.AtomsFM(parameters)


## Use turtle plotter for running plots ###
plot  = xyplotter.xyplotter(0.0,0.0,1.0,1.0, "Plot of Feature Functions",4)



##### read in atoms lazilly - For debugging only generate features functions of one atom ####
#### plotting the logarithm of the feature functions ####
iatom = 0  # change the atom
data = ''
count = 0
for (symbols,rions,fions,energy) in read_fei_urlfile(feifilename):
   nion = len(symbols)
   framefm = []
   #for ii in range(nion):
   #   framefm.append(afm(rions,ii))
   framefm.append(afm(rions,iatom))

   sstr = "%d " % count
   for i in range(52): sstr += "%f " %framefm[0][i]
   #for i in range(8): sstr += "%f " %framefm[0][i]
   print(sstr)
   sstrl = "%d " % count
   for i in range(52): sstrl += "%f " % (math.log(framefm[0][i]))
   data += sstrl + '\n'
   #print("frame number,framefm=",count,framefm[0][0:9],framefm[1])
   if ((count%50)==0): plot_ffpaths(plot,data)
   count += 1



### plot the relative energies using html ###
data = "#Title A simpler html plot for Raymond\n" + data
html = xydata_plotdatajs(data,"XY plot")
path = os.path.abspath('fei2ff.html')
url = 'file://' + path
with open(path, 'w') as ff:
   ff.write(html)
webbrowser.open(url)

###wait for return so that plot can be seen###
x = raw_input("--Press return to finish--")
print("yeh")
