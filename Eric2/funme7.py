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
parameterfilename = "../Eric/parameters0b"
parameters = []
with open(parameterfilename,'r') as ff: 
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
feidatafile = "../Eric/h2-aimd.fei"
with open(feidatafile,'r') as ff: 
   feidata = ff.read()
nion = int(feidata.split('\n')[0])
print("nion=",nion)
count = 0
feidata = feidata.split('\n')
nframes = len(feidata)/(nion+2+3)
print("nframes,nion,filelines = ",nframes,nion,len(feidata)-1)
energylist = []
rionslist  = []
fionslist  = []
#tobohr = 1.0/0.529177
tobohr = 1.0
for frame in range(nframes):
   energylist.append(float(feidata[(nion+5)*frame+1]))
   geom = feidata[(nion+5)*frame+5:(nion+5)*(frame+1)]
   rions  = [tobohr*float(ln.split()[i]) for ln in geom for i in range(3,6)]
   fions  = [tobohr*float(ln.split()[i]) for ln in geom for i in range(6,9)]
   rionslist.append(rions)
   fionslist.append(fions)


#### define Atomic Feature Mapping ####
afm = atomfm.AtomsFM(parameters)

#### define NN machine and weights for all atoms ####
nparameters = len(parameters)
sigmoid     = lambda x: 1.0/(1.0+math.exp(-x))
sigmoidp    = lambda x: math.exp(-x)/(math.exp(-x)+1.0)**2
sigmoidpp   = lambda x: 2*math.exp(-2*x)/(math.exp(-x)+1.0)**3 - math.exp(-x)/(math.exp(-x)+1.0)**2
nn_machine  = myfeedforward.MyFeedForward([nparameters,nparameters,1],[sigmoid,sigmoid,lambda x: x],[sigmoidp,sigmoidp,lambda x: 1], [sigmoidpp,sigmoidpp,lambda x: 0])
nn_weights  = []
for ii in range(nion):
   nn_weights += nn_machine.initial_w()


print("lenght w=",len(nn_weights))
nw = len(nn_weights)/nion
print("nw=",nw)

   
#nframes=200
for it0 in range(10):
   alpha = 0.1/nion
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

      #print("it0,frame,nion,len(framefm),eguess,exact,error=",it0,frame,nion,len(framefm),sum(etmp),energylist[frame],error,math.sqrt(error))
      print(" ",it0*nframes+frame,nion,len(framefm),sum(etmp),energylist[frame],error,math.sqrt(error))
      allframesfm.append(framefm)
      for i in range(len(nn_weights)):
         nn_weights[i] -= alpha*derror1detmp*dedw[i]

print("checking Energies and Forces")
nion3 = 3*nion
for frame in range(nframes):
   esum = 0.0
   force3 = [0.0]*nion3
   for ii in range(nion):
      #tafm = afm(rionslist[frame],ii)
      fafm = afm.Egradients(rionslist[frame],ii)
      eee = nn_machine.evaluate(fafm[0],nn_weights[ii*nw:(ii+1)*nw])
      esum += eee[0]
      fff = nn_machine.gradients_evaluate(fafm[0],nn_weights[ii*nw:(ii+1)*nw])
      #print("fff=",len(fff),nparameters,fff)

      for jj in range(nion3):
         for k in range(nparameters):
            force3[jj] -= fafm[1][jj + k*nion3]*fff[k]

   print("frame,esum=",frame,esum,energylist[frame],force3,fionslist[frame])


rion0 = [0.0]*nion3
rion1 = [0.0]*nion3
rion2 = [0.0]*nion3
vion  = [0.0]*nion3
for ii in range(nion3):
   rion0[ii] = rionslist[0][ii]
   rion1[ii] = rionslist[1][ii]
   rion2[ii] = rionslist[2][ii]

dt = 5.0
mass = 1822.89*2.0
dataxyz = open("funme7.xyz",'w')
datae   = open("funme7.emotion",'w')

for tt in range(1000):
   for ii in range(nion3):
      rion0[ii] = rion1[ii]
      rion1[ii] = rion2[ii]

   vsum   = 0.0
   force3 = [0.0]*nion3
   for ii in range(nion):
      fafm = afm.Egradients(rion1,ii)
      eee = nn_machine.evaluate(fafm[0],nn_weights[ii*nw:(ii+1)*nw])
      fff = nn_machine.gradients_evaluate(fafm[0],nn_weights[ii*nw:(ii+1)*nw])
      vsum += eee[0]
      for jj in range(nion3):
         for k in range(nparameters):
            force3[jj] -= fafm[1][jj + k*nion3]*fff[k]

   for ii in range(nion3):
      rion2[ii] = 2*rion1[ii] - rion0[ii] + (dt*dt/mass)*force3[ii]
      #print("ii,rion0,rion1,rion2,force3=",ii,rion0[ii],rion1[ii],rion2[ii],force3[ii])

   for ii in range(nion3):
      vion[ii] = (rion2[ii]-rion0[ii])/(2.0*dt)

   ke = 0.0
   for ii in range(nion3):
      ke += 0.5*mass*vion[ii]*vion[ii]

   print(tt*dt,ke+vsum,vsum,ke)
   datae.write("%f %f %f %f\n" % (tt*dt,ke+vsum,vsum,ke))

   ssframe = "%d\n\n" % (nion)
   for ii in range(nion):
      ssframe += "%s %f %f %f\n" % ('H',rion1[3*ii]*0.529177,rion1[3*ii+1]*0.529177,rion1[3*ii+2]*0.529177)
   dataxyz.write(ssframe)

dataxyz.close()
datae.close()
exit()

print("len(allframesfm)=",len(allframesfm))
print("\n\n")
datafp = open("fframeerror4a.dat",'w')


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
