import math

class AtomFunctions:
   """AtomFunctions - defines functions for Atom Feature mapping

   The following functions are defined
      1) f1(rion,i) = Sum(j!=i) fcut(rij)
      2) f2(rion,i) = Sum(j!=i) exp(-eta*(rij-rs)**2)*fcut(rij)
      3) f3(rion,i) = Sum(j!=i) cos(kappa*rij)*fcut(rij)
      4) f4(rion,i) = 2**(1-xi) * Sum(j!=i) Sum(k!=i,k!=j) (1-lmbda*cos(thetaijk))**xi exp(-eta*(rij*rij + rik*rik + rjk*rjk))*fcut(rij)*fcut(rik)*fcut(rjk)
      5) f5(rion,i) = 2**(1-xi) * Sum(j!=i) Sum(k!=i) (1-lmbda*cos(thetaijk))**xi exp(-eta*(rij*rij + rik*rik))*fcut(rij)*fcut(rik)

   Defining functions  
      f1 = AtomFunctions([1,rcut])
      f2 = AtomFunctions([2,rcut,eta,rs])
      f3 = AtomFunctions([3,rcut,kappa])
      f4 = AtomFunctions([4,rcut,eta,lmbda,xi])
      f5 = AtomFunctions([4,rcut,eta,lmbda,xi])

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

   def __f2exp(self,eta,rs,r):
      return math.exp(-eta*(r-rs)*(r-rs))

   def __f2cos(self,kappa,r):
      return math.cos(kappa*r)

   def __dfcutdr(self,rcut,r):
      ff = 0.0
      if r < rcut:
         ff = -0.5*math.pi/5.0*math.sin(math.pi*x/5.0)
      return ff

   def __df2expdr(self,eta,rs,r):
      return -2.0*eta*(r-rs)*math.exp(-eta*(r-rs)*(r-rs))

   def __df2cosdr(self,kappa,r):
      return -kappa*math.sin(kappa*r)


   def __dist(self,rion,i,j):
      x = rion[3*j]   - rion[3*i]
      y = rion[3*j+1] - rion[3*i+1]
      z = rion[3*j+2] - rion[3*i+2]
      return math.sqrt(x*x + y*y + z*z)

   def __grad_dist(self,rion,i,j):
      x = rion[3*j]   - rion[3*i]
      y = rion[3*j+1] - rion[3*i+1]
      z = rion[3*j+2] - rion[3*i+2]
      r = math.sqrt(x*x + y*y + z*z)
      fi = [-x/r,-y/r,-z/r]
      fj = [ x/r, y/r, z/r]
      return (fi,fj)

   def __theta(self,rion,i,j,k):
      if (j==k):
         vv = 0.0
      else:
         xij = rion[3*i]   - rion[3*j]
         yij = rion[3*i+1] - rion[3*j+1]
         zij = rion[3*i+2] - rion[3*j+2]
         xik = rion[3*i]   - rion[3*k]
         yik = rion[3*i+1] - rion[3*k+1]
         zik = rion[3*i+2] - rion[3*k+2]
         rij = math.sqrt(xij*xij+yij*yij+zij*zij)
         rik = math.sqrt(xik*xik+yik*yik+zik*zik)
         ff =(xij*xik + yij*yik + zij*zik)/(rij*rik)
         if (ff>1.0):    ff = 1.0
         if (ff<(-1.0)): ff = -1.0
         vv = math.acos(ff)
      return vv

   def __grad_theta(self,rion,i,j,k):
      xij = rion[3*i]   - rion[3*j]
      yij = rion[3*i+1] - rion[3*j+1]
      zij = rion[3*i+2] - rion[3*j+2]
      xik = rion[3*i]   - rion[3*k]
      yik = rion[3*i+1] - rion[3*k+1]
      zik = rion[3*i+2] - rion[3*k+2]
      rij = math.sqrt(xij*xij+yij*yij+zij*zij)
      rik = math.sqrt(xik*xik+yik*yik+zik*zik)
      ctheta =(xij*xik + yij*yik + zij*zik)/(rij*rik)
      if (ctheta>1.0):    ctheta = 1.0
      if (ctheta<(-1.0)): ctheta = -1.0
      stheta = math.sqrt(1.0-ctheta*ctheta)
      if (stheta<0.001): stheta=0.001
      stheta = 1.0/stheta
      q = math.acos(ctheta)
      fi = []
      fj = []
      fk = []
      return (fi,fj,fk)

   def __init__(self,parameters):
      self.__ftype = parameters[0]
      self.__parameters = parameters[1:]
      return

   def __call__(self,rion,i):
      nion = int(len(rion)/3)
      ff = 0.0
      if self.__ftype==1:
         rcut = self.__parameters[0]
         for j in range(nion):
            if (i!=j):
               r = self.__dist(rion,i,j)
               ff += self.__fcut(rcut,r)
      elif self.__ftype==2:
         rcut = self.__parameters[0]
         eta  = self.__parameters[1]
         rs   = self.__parameters[2]
         for j in range(nion):
            if (i!=j):
               r = self.__dist(rion,i,j)
               ff += self.__f2exp(eta,rs,r)*self.__fcut(rcut,r)
      elif self.__ftype==3:
         rcut  = self.__parameters[0]
         kappa = self.__parameters[1]
         for j in range(nion):
            if (i!=j):
               r = self.__dist(rion,i,j)
               ff += self.__f2cos(kappa,r)*self.__fcut(rcut,r)
      elif self.__ftype==4:
         rcut  = self.__parameters[0]
         eta   = self.__parameters[1]
         lmbda = self.__parameters[2]
         xi    = self.__parameters[3]
         for j in range(nion):
            if (i!=j):
               rij  = self.__dist(rion,i,j)
               rij2 = rij*rij
               for k in range(nion):
                  if (i!=k) and (j!=k):
                     rik     = self.__dist(rion,i,k)
                     rjk     = self.__dist(rion,j,k)
                     rik2 = rik*rik
                     rjk2 = rjk*rjk
                     theta = self.__theta(rion,i,j,k)
                     tmp1 = (1.0-lmbda*math.cos(theta))**xi
                     tmp2 = math.exp(-eta*(rij2 + rik2 + rjk2))
                     ff += tmp1*tmp2*self.__fcut(rcut,rij)*self.__fcut(rcut,rik)*self.__fcut(rcut,rjk)
         ff *= 2.0**(1-xi)
      elif self.__ftype==5:
         rcut  = self.__parameters[0]
         eta   = self.__parameters[1]
         lmbda = self.__parameters[2]
         xi    = self.__parameters[3]
         for j in range(nion):
            if (i!=j):
               rij  = self.__dist(rion,i,j)
               rij2 = rij*rij
               for k in range(nion):
                  if (i!=k) and (j!=k):
                     rik     = self.__dist(rion,i,k)
                     rjk     = self.__dist(rion,j,k)
                     rik2 = rik*rik
                     rjk2 = rjk*rjk
                     theta = self.__theta(rion,i,j,k)
                     tmp1 = (1.0-lmbda*math.cos(theta))**xi
                     tmp2 = math.exp(-eta*(rij2 + rik2))
                     ff += tmp1*tmp2*self.__fcut(rcut,rij)*self.__fcut(rcut,rik)
         ff *= 2.0**(1-xi)
      else:
         print("no function defined, ftype=",self.__ftype)
     
      return ff

class AtomsFM:

   def __fcut(self,rcut,r):
      ff = 0.0
      if r < rcut:
         ff = 0.5*(math.cos(math.pi*r/rcut)+1.0)
      return ff

   def __f2exp(self,eta,rs,r):
      return math.exp(-eta*(r-rs)*(r-rs))

   def __f2cos(self,kappa,r):
      return math.cos(kappa*r)

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
      fm = [0.0]*self.fmsize
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
               for p2 in self.p2s:
                  fm[p2[0]] += self.__f2exp(p2[2],p2[3],rij)*self.__fcut(p2[1],rij)*(0.5*p2[2]/math.pi)**0.5 
               for p3 in self.p3s:
                  fm[p3[0]] += self.__f2cos(p3[2],rij)*self.__fcut(p3[1],rij)

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
                        for p5 in self.p5s:
                           tmp0 = 2.0**(1-p5[4])
                           tmp1 = (1.0-p5[3]*ctheta)**p5[4]
                           tmp2 = math.exp(-p5[2]*(rij2 + rik2)) * (0.5*p5[2]/math.pi)**0.5 
                           fm[p5[0]] += tmp0*tmp1*tmp2*self.__fcut(p5[1],rij)*self.__fcut(p5[1],rik)
      return fm

   def gradients(self,rion,i):
      nion = int(len(rion)/3)
      gradfm = [0.0]*len(rion)*self.fmsize

      return gradrm


