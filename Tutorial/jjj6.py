import myfeedforward5 as myfeedforward
import math,random
import xyplotter

def plote(plot,A,machine,weights,weights0):

   x1 = []
   y  = []
   y0 = []
   y1 = []
   dy = 0.1
   ymin = 99999.9
   ymax = -99999.9
   for i in range(-19,20):
      x = i/(19.0)*2
      e = A*(x*x - 0.1*x**6)
      e1 = machine.evaluate([x],weights)[0]
      e0 = machine.evaluate([x],weights0)[0]

      if (e  < ymin): ymin = e
      if (e0 < ymin): ymin = e0
      if (e1 < ymin): ymin = e1
      if (e  > ymax): ymax = e
      if (e0 > ymax): ymax = e0
      if (e1 > ymax): ymax = e1

      x1.append(x)
      y.append(e)
      y0.append(e0)
      y1.append(e1)
   plot.resetwindow(-2.0,ymin-dy,2.0,ymax+dy,"Spring Energies")
   plot.plot(x1,y,"black")
   plot.plot(x1,y0,"red")
   plot.plot(x1,y1,"blue")

alpha0 = 0.01
alpha = 0.001
beta1 = 0.9
beta2 = 0.999
eps   = 1e-8


A = 0.2

beta = 1.0
#sigmoid   = lambda x: 1.0/(1.0+math.exp(-x))
#sigmoidp  = lambda x: math.exp(-x)/(1.0+math.exp(-x))**2
#sigmoidpp = lambda x: math.exp(-x)*(math.exp(-x)-1.0)/(1.0+math.exp(-x))**3
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
machine = myfeedforward.MyFeedForward([1,13,26,52,1],[xmoid1,sigmoid,sigmoid,sigmoid,xmoid1],[xmoidp1,sigmoidp,sigmoidp,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp1],bias)
#machine = myfeedforward.MyFeedForward([1,52,1],[xmoid1,sigmoid,xmoid1],[xmoidp1,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,xmoidpp1],bias)

#machine = myfeedforward.MyFeedForward([1,1,1,1,1,1,1],[sigmoid,sigmoid,sigmoid,sigmoid,xmoid3,xmoid2,xmoid1],[sigmoidp,sigmoidp,sigmoidp,sigmoidp,xmoidp3,xmoidp2,xmoidp1],[sigmoidpp,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp3,xmoidpp2,xmoidpp1],bias)

#machine = myfeedforward.MyFeedForward([1,1,1,1,1,1],[sigmoid,sigmoid,sigmoid,sigmoid,xmoid2,xmoid1],[sigmoidp,sigmoidp,sigmoidp,sigmoidp,xmoidp2,xmoidp1],[sigmoidpp,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp2,xmoidpp1],bias)
weights = machine.initial_w()

for i in range(len(weights)):
   weights[i] *= 1.0

weights0 = weights[:]
m = [0.0]*len(weights)
v = [0.0]*len(weights)
beta1t = 1.0
beta2t = 1.0

plot  = xyplotter.xyplotter(-1.0,0.0,1.0,1.0, "Spring Energies",4)

plote(plot,A,machine,weights,weights0)

enter0 = raw_input(" -- start simulation --- ")

for i in range(1000000):

   x = 4.0*random.random()-2.0
   #e = A*x*x
   e = A*(x*x - 0.1*x**6)
   e1 = machine.evaluate([x],weights)[0]
   gg = machine.w_energy_gradient([x],[e],weights)
   error = gg[0]
   g1 = gg[1]
   beta1t *= beta1
   beta2t *= beta2
   alphat = alpha*math.sqrt(1.0-beta2t)/(1.0-beta1t)
   for j in range(len(g1)):
      m[j] = beta1*m[j] + (1.0-beta1)*g1[j]
      v[j] = beta2*v[j] + (1.0-beta2)*g1[j]*g1[j]
      weights[j] -= alphat*m[j]/(math.sqrt(v[j]) + eps)

   e0 = machine.evaluate([x],weights0)[0]
   gg0 = machine.w_energy_gradient([x],[e],weights0)
   error0 = gg0[0]
   g0 = gg0[1]
   for j in range(len(g0)):
      weights0[j] -= alpha0*g0[j]

   #print "x,e,e1,error,w=",x,e,e1,error,weights
   print x,e,e1,(e1-e)**2,error,"  |  ",e,e0,error0

   if ((i%1000)==0):
      plote(plot,A,machine,weights,weights0)

   #if ((i%50000)==0):
   #   beta *= 0.9
   #   if (beta<0.1): beta=0.1
   #   sigmoid      = lambda x: 0.5*(math.tanh(beta*x)+1.0)
   #   sigmoidp     = lambda x: 0.5*beta*(1.0/math.cosh(beta*x))**2
   #   sigmoidpp    = lambda x: 0.5*(-2.0)*beta*beta*math.tanh(beta*x)*(1.0/math.sech(beta*x))**2
   #   machine = myfeedforward.MyFeedForward([1,52,1],[xmoid1,sigmoid,xmoid1],[xmoidp1,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,xmoidpp1],bias)
   #   #machine = myfeedforward.MyFeedForward([1,52,1],[xmoid1,sigmoid,xmoid1],[xmoidp1,sigmoidp,xmoidp1],[xmoidpp1,sigmoidpp,xmoidpp1],bias)

enter1 = raw_input(" -- finished --- ")



