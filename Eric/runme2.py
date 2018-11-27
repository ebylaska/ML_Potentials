from __future__ import print_function
import math,atomfm,myfeedforward3

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


#### define Atomic Feature Mapping ####
afm = atomfm.AtomsFM(parameters)

#### define NN machine and weights for all atoms ####
nparameters = len(parameters)
sigmoid     = lambda x: 1.0/(1.0+math.exp(-x))
sigmoidp    = lambda x: math.exp(-x)/(math.exp(-x)+1.0)**2
sigmoidpp   = lambda x: 2*math.exp(-2*x)/(math.exp(-x)+1.0)**3 - math.exp(-x)/(math.exp(-x)+1.0)**2
nn_machine  = myfeedforward3.MyFeedForward([nparameters,nparameters,nparameters,1],[sigmoid,sigmoid,sigmoid,lambda x: x],[sigmoidp,sigmoidp,sigmoidp,lambda x: 1], [sigmoidpp,sigmoidpp,sigmoidpp,lambda x: 0])
nn_weights  = []
for ii in range(nion):
   nn_weights += nn_machine.initial_w()


print("lenght w=",len(nn_weights))
nw = len(nn_weights)/nion
print("nw=",nw)

   
alpha = 0.001
allframesfm = []
for frame in range(nframes):
   framefm = []
   etmp = []
   dedw = []
   for ii in range(nion):
      framefm.append(afm(rionslist[frame],ii))
      #etmp += nn_machine.evaluate(framefm[ii],nn_weights[ii*nw:(ii+1)*nw])[0]
      eee = nn_machine.dyoutdw_gradient(framefm[ii],nn_weights[ii*nw:(ii+1)*nw])
      etmp += eee[0]
      dedw += eee[1]
      #print("ii,len(eee)=",ii,len(eee),len(eee[0]),len(eee[1]))
      #fff = nn_machine.gradients_evaluate(framefm[ii],nn_weights[ii])
      #print("ii,len(fff)=",ii,len(fff))
      #ggg = nn_machine.ddyoutdxdw_gradient(framefm[ii],nn_weights[ii])
      #print("ii,len(ggg)=",ii,len(ggg),len(ggg[0]),len(ggg[1]))

   error = (sum(etmp) - energylist[frame])**2
   derror1detmp    = 2.0*(sum(etmp) - energylist[frame])

   print("frame,nion,len(framefm),eguess,exact,error=",frame,nion,len(framefm),sum(etmp),energylist[frame],error)
   allframesfm.append(framefm)
   for i in range(len(nn_weights)):
      nn_weights[i] -= alpha*derror1detmp*dedw[i]

print("len(allframesfm)=",len(allframesfm))

