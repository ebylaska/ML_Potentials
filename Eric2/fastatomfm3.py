import math


class AtomsFM:
   """AtomFunctions - defines functions for Atom Feature mapping

   The following functions are defined
      1) f1(rion,i) = Sum(j!=i) fcut(rij)
      2) f2(rion,i) = Sum(j!=i) exp(-eta*(rij-rs)**2)*fcut(rij)
      3) f3(rion,i) = Sum(j!=i) cos(kappa*rij)*fcut(rij)
      4) f4(rion,i) = 2**(1-xi) * Sum(j!=i) Sum(k!=i,j) (1-lmbda*cos(thetaijk))**xi exp(-eta*(rij*rij + rik*rik + rjk*rjk))*fcut(rij)*fcut(rik)*fcut(rjk)
      5) f5(rion,i) = 2**(1-xi) * Sum(j!=i) Sum(k!=i,j) (1-lmbda*cos(thetaijk))**xi exp(-eta*(rij*rij + rik*rik))*fcut(rij)*fcut(rik)

   Defining functions  
      f1 = AtomFunctions([1,rcut])
      f2 = AtomFunctions([2,rcut,eta,rs])
      f3 = AtomFunctions([3,rcut,kappa])
      f4 = AtomFunctions([4,rcut,eta,lmbda,xi])
      f5 = AtomFunctions([5,rcut,eta,lmbda,xi])

   Using functions
      value1 = f1(rion,i)
      value2 = f2(rion,i)
      value3 = f3(rion,i)
      value4 = f4(rion,i)
      value5 = f5(rion,i)
 
   """


   def __fcut(self,rcut,r):
      ff = 0.0
      if r < rcut:
         ff = 0.5*(math.cos(math.pi*r/rcut)+1.0)
      return ff

   def __dfcut(self,rcut,r):
      dff = 0.0
      if r < rcut:
         #ff = 0.5*(math.cos(math.pi*r/rcut)+1.0)
         dff = -0.5*math.sin(math.pi*r/rcut) * math.pi/rcut
      return dff

   def __f2exp(self,eta,rs,r):
      return math.exp(-eta*(r-rs)*(r-rs))

   def __df2exp(self,eta,rs,r):
      #return math.exp(-eta*(r-rs)*(r-rs))
      return math.exp(-eta*(r-rs)*(r-rs))*(-2.0*eta*(r-rs))

   def __f2cos(self,kappa,r):
      return math.cos(kappa*r)

   def __df2cos(self,kappa,r):
      #return math.cos(kappa*r)
      return -math.sin(kappa*r)*kappa

   def __init__(self,parameterslist):
      #print("hello???")
      self.p1s = []
      self.p2s = []
      self.p3s = []
      self.p4s = []
      self.p5s = []
      self.fmsize = len(parameterslist)
      count = 0
      for parameters in parameterslist:
         if parameters[0] == 1: self.p1s.append([count]+parameters[1:])
         if parameters[0] == 2: self.p2s.append([count]+parameters[1:])
         if parameters[0] == 3: self.p3s.append([count]+parameters[1:])
         if parameters[0] == 4: self.p4s.append([count]+parameters[1:])
         if parameters[0] == 5: self.p5s.append([count]+parameters[1:])
         count += 1
      #print("afm init fmsize=",self.fmsize)
      return

   def __call__(self,rion,i):
      #print("fmsize=",self.fmsize)
      #print("len rion=",len(rion))
      fm  = [0.0]*self.fmsize
      cnt = [0.0]*self.fmsize
      nion = int(len(rion)/3)
      #print(len(rion), nion)
      Rcut = 6.0/0.529177
      for j in range(nion):
         if (i!=j):
            xij = rion[3*i]   - rion[3*j]
            yij = rion[3*i+1] - rion[3*j+1]
            zij = rion[3*i+2] - rion[3*j+2]
            rij2 = xij*xij + yij*yij + zij*zij
            rij  = math.sqrt(rij2)
            if (rij<=Rcut):
               for p1 in self.p1s:
                  fm[p1[0]] += self.__fcut(p1[1],rij)
                  cnt[p1[0]]+= 1.0
               for p2 in self.p2s:
                  fm[p2[0]] += self.__f2exp(p2[2],p2[3],rij)*self.__fcut(p2[1],rij)*(0.5*p2[2]/math.pi)**0.5 
                  cnt[p2[0]]+= 1.0
               for p3 in self.p3s:
                  fm[p3[0]] += self.__f2cos(p3[2],rij)*self.__fcut(p3[1],rij)
                  cnt[p3[0]]+= 1.0

               for k in range(nion):
                  if (i!=k) and (j!=k):
                     xik = rion[3*i]   - rion[3*k]
                     yik = rion[3*i+1] - rion[3*k+1]
                     zik = rion[3*i+2] - rion[3*k+2]
                     xjk = rion[3*j]   - rion[3*k]
                     yjk = rion[3*j+1] - rion[3*k+1]
                     zjk = rion[3*j+2] - rion[3*k+2]
                     rik2 = xik*xik + yik*yik + zik*zik
                     rjk2 = xjk*xjk + yjk*yjk + zjk*zjk
                     rik  = math.sqrt(rik2)
                     rjk  = math.sqrt(rjk2)
                     if (rik<Rcut) and (rjk<=Rcut):
                        ctheta =(xij*xik + yij*yik + zij*zik)/(rij*rik)
                        if (ctheta>1.0):    ctheta = 1.0
                        if (ctheta<(-1.0)): ctheta = -1.0
                        for p4 in self.p4s:
                           tmp0 = 2.0**(1-p4[4])
                           tmp1 = (1.0-p4[3]*ctheta)**p4[4]
                           tmp2 = math.exp(-p4[2]*(rij2 + rik2 + rjk2)) * (0.2*p4[2]/math.pi)**0.5
                           fm[p4[0]] += tmp0*tmp1*tmp2*self.__fcut(p4[1],rij)*self.__fcut(p4[1],rik)*self.__fcut(p4[1],rjk)
                           cnt[p4[0]]+= 1.0
                        for p5 in self.p5s:
                           tmp0 = 2.0**(1-p5[4])
                           tmp1 = (1.0-p5[3]*ctheta)**p5[4]
                           tmp2 = math.exp(-p5[2]*(rij2 + rik2)) * (0.5*p5[2]/math.pi)**0.5 
                           fm[p5[0]] += tmp0*tmp1*tmp2*self.__fcut(p5[1],rij)*self.__fcut(p5[1],rik)
                           cnt[p5[0]]+= 1.0

      for p1 in self.p1s:
         if (cnt[p1[0]]>0.0):
            fm[p1[0]] /= cnt[p1[0]]
      for p2 in self.p2s:
         if (cnt[p2[0]]>0.0):
            fm[p2[0]] /= cnt[p2[0]]
      for p3 in self.p3s:
         if (cnt[p3[0]]>0.0):
            fm[p3[0]] /= cnt[p3[0]]
      for p4 in self.p4s:
         if (cnt[p4[0]]>0.0):
            fm[p4[0]] /= cnt[p4[0]]
      for p5 in self.p5s:
         if (cnt[p5[0]]>0.0):
            fm[p5[0]] /= cnt[p5[0]]

      return fm

   def Egradients(self,rion,i):
      nion3= int(len(rion))
      nion = int(len(rion)/3)
      fm     = [0.0]*self.fmsize
      cnt    = [0.0]*self.fmsize
      gradfm = [0.0]*nion3*self.fmsize

      Rcut = 6.0/0.529177
      for j in range(nion):
         if (i!=j):
            xij = rion[3*i]   - rion[3*j]
            yij = rion[3*i+1] - rion[3*j+1]
            zij = rion[3*i+2] - rion[3*j+2]
            rij2 = xij*xij + yij*yij + zij*zij
            rij  = math.sqrt(rij2)
            if (rij<=Rcut):
               for p1 in self.p1s:
                  fm[p1[0]] += self.__fcut(p1[1],rij)
                  dfmdr = self.__dfcut(p1[1],rij)
                  gradfm[p1[0]*nion3 + 3*i]     += dfmdr*(xij/rij)
                  gradfm[p1[0]*nion3 + 3*j]     -= dfmdr*(xij/rij)
                  gradfm[p1[0]*nion3 + 3*i + 1] += dfmdr*(yij/rij)
                  gradfm[p1[0]*nion3 + 3*j + 1] -= dfmdr*(yij/rij)
                  gradfm[p1[0]*nion3 + 3*i + 2] += dfmdr*(zij/rij)
                  gradfm[p1[0]*nion3 + 3*j + 2] -= dfmdr*(zij/rij)
                  cnt[p1[0]]+= 1.0
               for p2 in self.p2s:
                  fm[p2[0]] += self.__f2exp(p2[2],p2[3],rij)*self.__fcut(p2[1],rij)*(0.5*p2[2]/math.pi)**0.5
                  dfmdr = (self.__df2exp(p2[2],p2[3],rij)*self.__fcut(p2[1],rij) + self.__f2exp(p2[2],p2[3],rij)*self.__dfcut(p2[1],rij))*(0.5*p2[2]/math.pi)**0.5
                  gradfm[p2[0]*nion3 + 3*i]     += dfmdr*(xij/rij)
                  gradfm[p2[0]*nion3 + 3*j]     -= dfmdr*(xij/rij)
                  gradfm[p2[0]*nion3 + 3*i + 1] += dfmdr*(yij/rij)
                  gradfm[p2[0]*nion3 + 3*j + 1] -= dfmdr*(yij/rij)
                  gradfm[p2[0]*nion3 + 3*i + 2] += dfmdr*(zij/rij)
                  gradfm[p2[0]*nion3 + 3*j + 2] -= dfmdr*(zij/rij)
                  cnt[p2[0]]+= 1.0
               for p3 in self.p3s:
                  fm[p3[0]] += self.__f2cos(p3[2],rij)*self.__fcut(p3[1],rij)
                  dfmdr = (self.__df2cos(p3[2],rij)*self.__fcut(p3[1],rij) + self.__f2cos(p3[2],rij)*self.__dfcut(p3[1],rij))
                  gradfm[p3[0]*nion3 + 3*i]     += dfmdr*(xij/rij)
                  gradfm[p3[0]*nion3 + 3*j]     -= dfmdr*(xij/rij)
                  gradfm[p3[0]*nion3 + 3*i + 1] += dfmdr*(yij/rij)
                  gradfm[p3[0]*nion3 + 3*j + 1] -= dfmdr*(yij/rij)
                  gradfm[p3[0]*nion3 + 3*i + 2] += dfmdr*(zij/rij)
                  gradfm[p3[0]*nion3 + 3*j + 2] -= dfmdr*(zij/rij)
                  cnt[p3[0]]+= 1.0

               for k in range(nion):
                  if (i!=k) and (j!=k):
                     xik = rion[3*i]   - rion[3*k]
                     yik = rion[3*i+1] - rion[3*k+1]
                     zik = rion[3*i+2] - rion[3*k+2]
                     xjk = rion[3*j]   - rion[3*k]
                     yjk = rion[3*j+1] - rion[3*k+1]
                     zjk = rion[3*j+2] - rion[3*k+2]
                     rik2 = xik*xik + yik*yik + zik*zik
                     rjk2 = xjk*xjk + yjk*yjk + zjk*zjk
                     rik  = math.sqrt(rik2)
                     rjk  = math.sqrt(rjk2)
                     if (rik<Rcut) and (rjk<=Rcut):
                        ctheta =(xij*xik + yij*yik + zij*zik)/(rij*rik)
                        if (ctheta>1.0):    ctheta = 1.0
                        if (ctheta<(-1.0)): ctheta = -1.0
                        dcthetadxi = (xij+xik)/(rij*rik) - ctheta*(xij/rij**2 + xik/rik**2)
                        dcthetadyi = (yij+yik)/(rij*rik) - ctheta*(yij/rij**2 + yik/rik**2)
                        dcthetadzi = (zij+zik)/(rij*rik) - ctheta*(zij/rij**2 + zik/rik**2)
                        dcthetadxj = -xik/(rij*rik) + ctheta*(xij/rij**2)
                        dcthetadyj = -yik/(rij*rik) + ctheta*(yij/rij**2)
                        dcthetadzj = -zik/(rij*rik) + ctheta*(zij/rij**2)
                        dcthetadxk = -xij/(rij*rik) + ctheta*(xik/rik**2)
                        dcthetadyk = -yij/(rij*rik) + ctheta*(yik/rik**2)
                        dcthetadzk = -zij/(rij*rik) + ctheta*(zik/rik**2)

                        for p4 in self.p4s:
                           tmp0 = 2.0**(1-p4[4])
                           tmp1 = (1.0-p4[3]*ctheta)**p4[4]
                           tmp1b= (1.0-p4[3]*ctheta)**(p4[4]-1)
                           tmp2 = math.exp(-p4[2]*(rij2 + rik2 + rjk2)) * (0.2*p4[2]/math.pi)**0.5
                           fij = self.__fcut(p4[1],rij)
                           fik = self.__fcut(p4[1],rik)
                           fjk = self.__fcut(p4[1],rjk)
                           #ff0  = tmp0*tmp1*tmp2*self.__fcut(p4[1],rij)*self.__fcut(p4[1],rik)*self.__fcut(p4[1],rjk)
                           fm[p4[0]] += tmp0*tmp1*tmp2*fij*fik*fjk
                           
                           dfij = self.__dfcut(p4[1],rij)
                           dfik = self.__dfcut(p4[1],rik)
                           dfjk = self.__dfcut(p4[1],rjk)

                           gradfm[p4[0]*nion3 + 3*i]   += tmp0*(-p4[4]*p4[3]*tmp1b*dcthetadxi - tmp1*2.0*p4[2]*(xij + xik))*tmp2*fij*fik*fjk
                           gradfm[p4[0]*nion3 + 3*i]   += tmp0*tmp1*tmp2*(fik*dfij*xij/rij + fij*dfik*xik/rik)*fjk
                           gradfm[p4[0]*nion3 + 3*i+1] += tmp0*(-p4[4]*p4[3]*tmp1b*dcthetadyi - tmp1*2.0*p4[2]*(yij + yik))*tmp2*fij*fik*fjk
                           gradfm[p4[0]*nion3 + 3*i+1] += tmp0*tmp1*tmp2*(fik*dfij*yij/rij + fij*dfik*yik/rik)*fjk
                           gradfm[p4[0]*nion3 + 3*i+2] += tmp0*(-p4[4]*p4[3]*tmp1b*dcthetadzi - tmp1*2.0*p4[2]*(zij + zik))*tmp2*fij*fik*fjk
                           gradfm[p4[0]*nion3 + 3*i+2] += tmp0*tmp1*tmp2*(fik*dfij*zij/rij + fij*dfik*zik/rik)*fjk

                           gradfm[p4[0]*nion3 + 3*j]   += tmp0*(-p4[4]*p4[3]*tmp1b*dcthetadxj + tmp1*2.0*p4[2]*(xij - xjk))*tmp2*fij*fik*fjk
                           gradfm[p4[0]*nion3 + 3*j]   += tmp0*tmp1*tmp2*(-fjk*dfij*xij/rij + fij*dfjk*xjk/rjk)*fik
                           gradfm[p4[0]*nion3 + 3*j+1] += tmp0*(-p4[4]*p4[3]*tmp1b*dcthetadyj + tmp1*2.0*p4[2]*(yij - yjk))*tmp2*fij*fik*fjk
                           gradfm[p4[0]*nion3 + 3*j+1] += tmp0*tmp1*tmp2*(-fjk*dfij*yij/rij + fij*dfjk*yjk/rjk)*fik
                           gradfm[p4[0]*nion3 + 3*j+2] += tmp0*(-p4[4]*p4[3]*tmp1b*dcthetadzj + tmp1*2.0*p4[2]*(zij - zjk))*tmp2*fij*fik*fjk
                           gradfm[p4[0]*nion3 + 3*j+2] += tmp0*tmp1*tmp2*(-fjk*dfij*zij/rij + fij*dfjk*zjk/rjk)*fik

                           gradfm[p4[0]*nion3 + 3*k]   += tmp0*(-p4[4]*p4[3]*tmp1b*dcthetadxk + tmp1*2.0*p4[2]*(xik + xjk))*tmp2*fij*fik*fjk
                           gradfm[p4[0]*nion3 + 3*k]   += tmp0*tmp1*tmp2*(-fjk*dfik*xik/rik - fik*dfjk*xjk/rjk)*fij
                           gradfm[p4[0]*nion3 + 3*k+1] += tmp0*(-p4[4]*p4[3]*tmp1b*dcthetadyk + tmp1*2.0*p4[2]*(yik + yjk))*tmp2*fij*fik*fjk
                           gradfm[p4[0]*nion3 + 3*k+1] += tmp0*tmp1*tmp2*(-fjk*dfik*yik/rik - fik*dfjk*yjk/rjk)*fij
                           gradfm[p4[0]*nion3 + 3*k+2] += tmp0*(-p4[4]*p4[3]*tmp1b*dcthetadzk + tmp1*2.0*p4[2]*(zik + zjk))*tmp2*fij*fik*fjk
                           gradfm[p4[0]*nion3 + 3*k+2] += tmp0*tmp1*tmp2*(-fjk*dfik*zik/rik - fik*dfjk*zjk/rjk)*fij

                           cnt[p4[0]]+= 1.0
                        for p5 in self.p5s:
                           tmp0 = 2.0**(1-p5[4])
                           tmp1 = (1.0-p5[3]*ctheta)**p5[4]
                           tmp1b= (1.0-p5[3]*ctheta)**(p5[4]-1)
                           tmp2 = math.exp(-p5[2]*(rij2 + rik2)) * (0.5*p5[2]/math.pi)**0.5
                           fij = self.__fcut(p5[1],rij)
                           fik = self.__fcut(p5[1],rik)
                           #fm[p5[0]] += tmp0*tmp1*tmp2*self.__fcut(p5[1],rij)*self.__fcut(p5[1],rik)
                           fm[p5[0]] += tmp0*tmp1*tmp2*fij*fik

                           dfij = self.__dfcut(p5[1],rij)
                           dfik = self.__dfcut(p5[1],rik)

                           gradfm[p5[0]*nion3 + 3*i]   += tmp0*(-p5[4]*p5[3]*tmp1b*dcthetadxi - tmp1*2.0*p5[2]*(xij + xik))*tmp2*fij*fik
                           gradfm[p5[0]*nion3 + 3*i]   += tmp0*tmp1*tmp2*(fik*dfij*xij/rij + fij*dfik*xik/rik)
                           gradfm[p5[0]*nion3 + 3*i+1] += tmp0*(-p5[4]*p5[3]*tmp1b*dcthetadyi - tmp1*2.0*p5[2]*(yij + yik))*tmp2*fij*fik
                           gradfm[p5[0]*nion3 + 3*i+1] += tmp0*tmp1*tmp2*(fik*dfij*yij/rij + fij*dfik*yik/rik)
                           gradfm[p5[0]*nion3 + 3*i+2] += tmp0*(-p5[4]*p5[3]*tmp1b*dcthetadzi - tmp1*2.0*p5[2]*(zij + zik))*tmp2*fij*fik
                           gradfm[p5[0]*nion3 + 3*i+2] += tmp0*tmp1*tmp2*(fik*dfij*zij/rij + fij*dfik*zik/rik)

                           gradfm[p5[0]*nion3 + 3*j]   += tmp0*(-p5[4]*p5[3]*tmp1b*dcthetadxj + tmp1*2.0*p5[2]*(xij))*tmp2*fij*fik
                           gradfm[p5[0]*nion3 + 3*j]   += tmp0*tmp1*tmp2*(-fik*dfij*xij/rij)
                           gradfm[p5[0]*nion3 + 3*j+1] += tmp0*(-p5[4]*p5[3]*tmp1b*dcthetadyj + tmp1*2.0*p5[2]*(yij))*tmp2*fij*fik
                           gradfm[p5[0]*nion3 + 3*j+1] += tmp0*tmp1*tmp2*(-fik*dfij*yij/rij)
                           gradfm[p5[0]*nion3 + 3*j+2] += tmp0*(-p5[4]*p5[3]*tmp1b*dcthetadzj + tmp1*2.0*p5[2]*(zij))*tmp2*fij*fik
                           gradfm[p5[0]*nion3 + 3*j+2] += tmp0*tmp1*tmp2*(-fik*dfij*zij/rij)

                           gradfm[p5[0]*nion3 + 3*k]   += tmp0*(-p5[4]*p5[3]*tmp1b*dcthetadxk + tmp1*2.0*p5[2]*(xik))*tmp2*fij*fik
                           gradfm[p5[0]*nion3 + 3*k]   += tmp0*tmp1*tmp2*(-fij*dfik*xik/rik)
                           gradfm[p5[0]*nion3 + 3*k+1] += tmp0*(-p5[4]*p5[3]*tmp1b*dcthetadyk + tmp1*2.0*p5[2]*(yik))*tmp2*fij*fik
                           gradfm[p5[0]*nion3 + 3*k+1] += tmp0*tmp1*tmp2*(-fij*dfik*yik/rik)
                           gradfm[p5[0]*nion3 + 3*k+2] += tmp0*(-p5[4]*p5[3]*tmp1b*dcthetadzk + tmp1*2.0*p5[2]*(zik))*tmp2*fij*fik
                           gradfm[p5[0]*nion3 + 3*k+2] += tmp0*tmp1*tmp2*(-fij*dfik*zik/rik)

                           cnt[p5[0]]+= 1.0

      for p1 in self.p1s:
         if (cnt[p1[0]]>0.0):
            fm[p1[0]] /= cnt[p1[0]]
            for jj in range(nion3): 
               gradfm[jj + p1[0]*nion3] /= cnt[p1[0]]
      for p2 in self.p2s:
         if (cnt[p2[0]]>0.0):
            fm[p2[0]] /= cnt[p2[0]]
            for jj in range(nion3): 
               gradfm[jj + p2[0]*nion3] /= cnt[p2[0]]
      for p3 in self.p3s:
         if (cnt[p3[0]]>0.0):
            fm[p3[0]] /= cnt[p3[0]]
            for jj in range(nion3): 
               gradfm[jj + p3[0]*nion3] /= cnt[p3[0]]
      for p4 in self.p4s:
         if (cnt[p4[0]]>0.0):
            fm[p4[0]] /= cnt[p4[0]]
            for jj in range(nion3): 
               gradfm[jj + p4[0]*nion3] /= cnt[p4[0]]
      for p5 in self.p5s:
         if (cnt[p5[0]]>0.0):
            fm[p5[0]] /= cnt[p5[0]]
            for jj in range(nion3): 
               gradfm[jj + p5[0]*nion3] /= cnt[p5[0]]

      return (fm,gradfm)


