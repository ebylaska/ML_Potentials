import myfeedforward5 as myfeedforward
import math,random
import xyplotter

def plote(plot,xdat,edat,machine,weights):

   x1 = []
   y  = []
   y1 = []
   dy = 0.1
   xmin = 99999.9
   xmax = -99999.9
   ymin = 99999.9
   ymax = -99999.9
   for i in range(len(xdat)):
      x = xdat[i]
      e = edat[i]
      e1 = machine.evaluate([x],weights)[0]

      if (x < xmin):  xmin = x
      if (x > xmax):  xmax = x
      if (e  < ymin): ymin = e
      if (e1 < ymin): ymin = e1
      if (e  > ymax): ymax = e
      if (e1 > ymax): ymax = e1

      x1.append(x)
      y.append(e)
      y1.append(e1)
   plot.resetwindow(xmin,ymin-dy,xmax,ymax+dy,"Spring Energies")
   plot.plot(x1,y,"black")
   plot.plot(x1,y1,"blue")

alpha0 = 0.01
alpha = 0.001
beta1 = 0.9
beta2 = 0.999
eps   = 1e-8


A = 0.2
xdat = []
edat = []
for i in range(50):
   x = 8.0*i/49.0
   e = A*(x*x - 0.1*x**6 + 0.5*(x-4)**8)
   xdat.append(x)
   edat.append(e)

beta = 1.0
#sigmoid   = lambda x: 1.0/(1.0+math.exp(-x))
#sigmoidp  = lambda x: math.exp(-x)/(1.0+math.exp(-x))**2
#sigmoidpp = lambda x: math.exp(-x)*(math.exp(-x)-1.0)/(1.0+math.exp(-x))**3
ap = 1.0
xp = 2.5
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
machine = myfeedforward.MyFeedForward([1,15,30,1],[xmoid1,sigmoid,sigmoid,xmoid1],[xmoidp1,sigmoidp,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,sigmoidpp,xmoidpp1],bias)
#machine = myfeedforward.MyFeedForward([1,13,26,52,1],[xmoid1,sigmoid,sigmoid,sigmoid,xmoid1],[xmoidp1,sigmoidp,sigmoidp,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp1],bias)
#machine = myfeedforward.MyFeedForward([1,52,1],[xmoid1,sigmoid,xmoid1],[xmoidp1,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,xmoidpp1],bias)

#machine = myfeedforward.MyFeedForward([1,1,1,1,1,1,1],[sigmoid,sigmoid,sigmoid,sigmoid,xmoid3,xmoid2,xmoid1],[sigmoidp,sigmoidp,sigmoidp,sigmoidp,xmoidp3,xmoidp2,xmoidp1],[sigmoidpp,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp3,xmoidpp2,xmoidpp1],bias)

#machine = myfeedforward.MyFeedForward([1,1,1,1,1,1],[sigmoid,sigmoid,sigmoid,sigmoid,xmoid2,xmoid1],[sigmoidp,sigmoidp,sigmoidp,sigmoidp,xmoidp2,xmoidp1],[sigmoidpp,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp2,xmoidpp1],bias)
weights = machine.initial_w()

for i in range(len(weights)):
   weights[i] *= 1.0

nw = len(weights)

print "weights=",weights
m = [0.0]*len(weights)
v = [0.0]*len(weights)
beta1t = 1.0
beta2t = 1.0

plot  = xyplotter.xyplotter(-1.0,0.0,1.0,1.0, "Spring Energies",4)

plote(plot,xdat,edat,machine,weights)

enter0 = raw_input(" -- start simulation --- ")

ii = 0
for i in range(1000000):

   #x = 4.0*random.random()-2.0
   #ii = random.randint(0,len(xdat)-1)
   x = xdat[ii]
   e = edat[ii]
   ii = ((ii+1) % len(xdat))

   gg = machine.w_energy_gradient([x],[e],weights)
   error = gg[0]
   g1 = gg[1]
   
   for j in range(1,nw-1):
      error += penalty(weights[j])
      g1[j] += penaltyp(weights[j])

   beta1t *= beta1
   beta2t *= beta2
   alphat = alpha*math.sqrt(1.0-beta2t)/(1.0-beta1t)
   for j in range(nw):
      m[j] = beta1*m[j] + (1.0-beta1)*g1[j]
      v[j] = beta2*v[j] + (1.0-beta2)*g1[j]*g1[j]
      weights[j] -= alphat*m[j]/(math.sqrt(v[j]) + eps)

   if ((i%10000)==0):
      e1 = machine.evaluate([x],weights)[0]
      print i,x,e,e1,error
      plote(plot,xdat,edat,machine,weights)


enter1 = raw_input(" -- finished --- ")



