from __future__ import print_function
import atomfm

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

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



#### read in atoms ####
with open("./eatomstructure.xyz",'r') as ff: 
   xyzdata = ff.read()
nion = int(xyzdata.split('\n')[0])
print("nion=",nion)
count = 0
xyzdata = xyzdata.split('\n')
nframes = len(xyzdata)/(nion+2)
print("nframes,nion,filelines = ",nframes,nion,len(xyzdata)-1)
energylist = []
rionslist  = []
fionslist  = []
tobohr = 1.0/0.529177
for frame in range(nframes):
   energylist.append(float(xyzdata[(nion+2)*frame+1]))
   geom = xyzdata[(nion+2)*frame+2:(nion+2)*(frame+1)]
   rions  = [tobohr*float(ln.split()[i]) for ln in geom for i in range(1,4)]
   fions  = [tobohr*float(ln.split()[i]) for ln in geom for i in range(4,7)]
   rionslist.append(rions)
   fionslist.append(fions)

#print("energy=",energylist,len(energylist))
#print("rionslist=",rionslist,len(rionslist))
#print("fionslist=",fionslist,len(fionslist))

   
for frame in range(nframes):
   framefm = []
   for ii in range(nion):
      framefm.append(afm(rionslist[frame],ii))

   print("frame,framefm=",frame,framefm)

