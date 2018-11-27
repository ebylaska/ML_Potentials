from __future__ import print_function
import math
import fastatomfm3 as atomfm
import myfeedforward4 as myfeedforward

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

#### read in parameters ####
parameters = []
with open("../Eric/parameters",'r') as ff: 
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
with open("../Eric/eatomstructure.xyz",'r') as ff: 
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
#tobohr = 1.0/0.529177
tobohr = 1.0
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
nn_machine  = myfeedforward.MyFeedForward([nparameters,nparameters,nparameters,1],[sigmoid,sigmoid,sigmoid,lambda x: x],[sigmoidp,sigmoidp,sigmoidp,lambda x: 1], [sigmoidpp,sigmoidpp,sigmoidpp,lambda x: 0])
nn_weights  = []
for ii in range(nion):
   nn_weights += nn_machine.initial_w()


print("lenght w=",len(nn_weights))
nw = len(nn_weights)/nion
print("nw=",nw)

   
nframes=200
alpha = 0.01/nion
allframesfm = []
for frame in range(nframes):
   framefm = []
   etmp = []
   dedw = []
   for ii in range(nion):
      framefm.append(afm(rionslist[frame],ii))
      #print("ldjfld?? =",frame,ii,rionslist[frame][3*ii],rionslist[frame][3*ii+1],rionslist[frame][3*ii+2])
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
print("\n\n")
datafp = open("frameerror3.dat",'w')


alpha /= nframes
for it in range(100):
   dedwall = [0.0]*len(nn_weights)
   errorall = 0.0
   errorall2 = 0.0
   if ((it%10)==0): print("     frame  error\n"+"      ----------------------")
   for frame in range(nframes):
      framefm = allframesfm[frame]
      etmp = []
      dedw = []
      for ii in range(nion):
         eee = nn_machine.dyoutdw_gradient(framefm[ii],nn_weights[ii*nw:(ii+1)*nw])
         etmp += eee[0]
         dedw += eee[1]
         #print("frame,ii,eee[0]=",frame,ii,eee[0])

      error = (sum(etmp) - energylist[frame])**2
      derrordetmp  = 2.0*(sum(etmp) - energylist[frame])
      if ((it%10)==0): 
         print("     ",frame,error,sum(etmp),energylist[frame])
      datafp.write("%d %10.6e %14.6f %14.6f\n" % (frame,error,sum(etmp),energylist[frame]))

      errorall += error
      errorall2 += error*error
      for i in range(len(nn_weights)):
         dedwall[i] += derrordetmp*dedw[i]

   if ((it%10)==0):
      print("      ----------------------")
   datafp.write("\n"); datafp.flush()
   ave1 = errorall/nframes
   ave2 = errorall2/nframes
   variance = ave2-ave1*ave1
   gerror = 0.0
   for i in range(len(nn_weights)):
      nn_weights[i] -= alpha*dedwall[i]
      gerror += dedwall[i]**2

   print("it, alpha, Average Error, Variance Error,gerror=",it,alpha,ave1,variance,math.sqrt(variance),gerror)
   print()


datafp.close()
