import myfeedforward5 as myfeedforward
import math,random
import xyplotter

def plote(plot,machine,weights):

   x0 = []
   z0 = []
   z1 = []
   z2 = []
   zz0 = []
   zz1 = []
   zz2 = []
   de = 0.1
   emin = 99999.99
   emax = -9999.99
   y0 = -1.0
   y1 =  0.5
   y2 =  0.0
   for i in range(-19,20):
      x = i/19.0
      e0 = math.sin(x*y0) + math.cos(x*y0)**2 + math.sin(x+y0)
      e1 = math.sin(x*y1) + math.cos(x*y1)**2 + math.sin(x+y1)
      e2 = math.sin(x*y2) + math.cos(x*y2)**2 + math.sin(x+y2)
      x0.append(x)
      z0.append(e0)
      z1.append(e1)
      z2.append(e2)
      ee0 = machine.evaluate([x,y0],weights)[0]
      ee1 = machine.evaluate([x,y1],weights)[0]
      ee2 = machine.evaluate([x,y2],weights)[0]
      zz0.append(ee0)
      zz1.append(ee1)
      zz2.append(ee2)
      if (e0 <emin): emin = e0
      if (ee0<emin): emin = ee0
      if (e1 <emin): emin = e1
      if (ee1<emin): emin = ee1
      if (e2 <emin): emin = e2
      if (ee2<emin): emin = ee2
      if (e0 >emax): emax = e0
      if (ee0>emax): emax = ee0
      if (e1 >emax): emax = e1
      if (ee1>emax): emax = ee1
      if (e2 >emax): emax = e2
      if (ee2>emax): emax = ee2

   plot.resetwindow(-1.0,emin-de,1.0,emax+de,"Spring Energies,xp=4.5")
   plot.plot(x0,z0, "blue")
   plot.dotplot(x0,zz0,"blue")

   plot.plot(x0,z1, "red")
   plot.dotplot(x0,zz1,"red")

   plot.plot(x0,z2, "green")
   plot.dotplot(x0,zz2,"green")



alpha0 = 0.01
alpha = 0.0001
beta1 = 0.9
beta2 = 0.999
eps   = 1e-8

beta = 2.0
#sigmoid   = lambda x: 1.0/(1.0+math.exp(-x))
#sigmoidp  = lambda x: math.exp(-x)/(1.0+math.exp(-x))**2
#sigmoidpp = lambda x: math.exp(-x)*(math.exp(-x)-1.0)/(1.0+math.exp(-x))**3
ap = 1.0
xp = 4.5
bp = 3.0
penalty  = lambda x: ap*(0.5*(math.tanh(bp*(x-xp)) - math.tanh(bp*(x+xp))) + 1.0)
penaltyp = lambda x: ap*0.5*bp*( (1/math.cosh(bp*(x-xp)))**2 - (1.0/math.cosh(bp*(x+xp)))**2)

sigmoid      = lambda x: 0.5*(math.tanh(beta*x)+1.0)
sigmoidp     = lambda x: 0.5*beta*(1.0/math.cosh(beta*x))**2
sigmoidpp    = lambda x: 0.5*(-2.0)*beta*beta*math.tanh(beta*x)*(1.0/math.sech(beta*x))**2
xmoid1   = lambda x: x
xmoidp1  = lambda x: 1.0
xmoidpp1 = lambda x: 0.0
#xmoid2   = lambda x: x*x * 0.5
#xmoidp2  = lambda x: x
#xmoidpp2 = lambda x: 1.0
#xmoid3   = lambda x: x*x*x * (1.0/6.0) 
#xmoidp3  = lambda x: 3*x*x * (1.0/6.0)
#xmoidpp3 = lambda x: 6*x   * (1.0/6.0)
#
#xmoid4   = lambda x: x*x*x*x * (1.0/24.0)
#xmoidp4  = lambda x: 4*x*x*x * (1.0/24.0)
#xmoidpp4 = lambda x: 12*x*x *  (1/0/24.0)

#bias = [[0.01],[0.01],[0.001],[0.0001],[0.00001],[0.0000001]]
bias = [[0.01],[0.01],[0.001]]
bias = []
machine = myfeedforward.MyFeedForward([2,20,40,1],[xmoid1,sigmoid,sigmoid,xmoid1],[xmoidp1,sigmoidp,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,sigmoidpp,xmoidpp1],bias)

weights = machine.initial_w()

for i in range(len(weights)):
   weights[i] *= 1.0

nw = len(weights)

print "weights=",weights
m = [0.0]*len(weights)
v = [0.0]*len(weights)
beta1t = 1.0
beta2t = 1.0

plot  = xyplotter.xyplotter(-1.0,0.0,1.0,1.0, "Spring Energies,xp=4.5",4)

plote(plot,machine,weights)

enter0 = raw_input(" -- start simulation --- ")

error = 999999.9
ii = 0
for i in range(1000000):

   xs = 2.0*random.random()-1.0
   ys = 2.0*random.random()-1.0
   es = math.sin(xs*ys) + math.cos(xs*ys)**2 + math.sin(xs+ys)

   gg = machine.w_energy_gradient([xs,ys],[es],weights)
   error0 = error
   error = gg[0]
   g1    = gg[1]

   for j in range(nw):
      error += penalty(weights[j])
      g1[j] += penaltyp(weights[j])


   beta1t *= beta1
   beta2t *= beta2
   alphat = alpha*math.sqrt(1.0-beta2t)/(1.0-beta1t)
   for j in range(nw):
      m[j] = beta1*m[j] + (1.0-beta1)*g1[j]
      v[j] = beta2*v[j] + (1.0-beta2)*g1[j]*g1[j]
      weights[j] -= alphat*m[j]/(math.sqrt(v[j]) + eps)

   if ((i%1000)==0):
      es1 = machine.evaluate([xs,ys],weights)[0]
      print "%10d %10.3f %10.3f %10.3f %10.3f  ||  %12.6f" % (i,xs,ys,es,es1,error)
      plote(plot,machine,weights)


enter1 = raw_input(" -- finished --- ")



