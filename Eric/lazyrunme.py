from __future__ import print_function
import atomfm

#xyzfilename = "small.xyz"
xyzfilename = "../eatomstructure.xyz"

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def read_in_lines(filename,nlines):
   """Lazy function (generator) to read a file in 
      nline chuncks."""
   with open(filename,'r') as fp:
      while True:
         data = ''.join([fp.readline() for _ in range(nlines)])
         if not data: 
            break
         yield data

def read_xyzfile(filename):
   """Lazy function (generator) to read a xyzfile 
      frame by frame."""
   tobohr = 1.0/0.529177
   with open(filename,'r') as ff: 
      xyzdata = ff.readline()
   nion = int(xyzdata.split('\n')[0])
   nlines = nion+2

   with open(filename,'r') as fp:
      while True:
         data = ''.join([fp.readline() for _ in range(nlines)])
         if not data: 
            break
         lines = data.strip().split('\n')
         energy  = float(lines[1])
         symbols = [ln.split()[0] for ln in lines[2:]]
         rions   = [tobohr*float(ln.split()[i]) for ln in lines[2:] for i in range(1,4)]
         fions   = [tobohr*float(ln.split()[i]) for ln in lines[2:] for i in range(4,7)]

         yield (symbols,rions,fions,energy)

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

print("parameters=",parameters,len(parameters))

#### define Atomic Feature Mapping ####
afm = atomfm.AtomsFM(parameters)



##### read in atoms lazilly ####
#with open("../eatomstructure.xyz",'r') as ff: 
#   xyzdata = ff.readline()
#nion = int(xyzdata.split('\n')[0])
#print("nion=",nion)

#count = 0
#tobohr = 1.0/0.529177
#for frame in read_in_lines("../eatomstructure.xyz",nion+2):
#   lines = frame.strip().split('\n')
#   energy = float(lines[1])
#   rions  = [tobohr*float(ln.split()[i]) for ln in lines[2:] for i in range(1,4)]
#   fions  = [tobohr*float(ln.split()[i]) for ln in lines[2:] for i in range(4,7)]
#
#   framefm = []
#   for ii in range(nion):
#      framefm.append(afm(rions,ii))
#
#   print("\n")
#   print("frame number,framefm=",count,framefm)
#   print("\n")
#   print("count=",count)
#   print("energy=",energy)
#   print("rions= ",rions)
#   print("fions= ",fions)
#   count += 1

count = 0
for (symbols,rions,fions,energy) in read_xyzfile(xyzfilename):
   nion = len(symbols)
   framefm = []
   for ii in range(nion):
      framefm.append(afm(rions,ii))

   print("\n")
   print("frame number,framefm=",count,framefm)
   print("\n")
   print("count=",count)
   print("nion=  ",nion)
   print("energy=",energy)
   print("rions= ",rions)
   print("fions= ",fions)
   count += 1

