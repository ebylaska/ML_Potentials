import myfeedforward5 as myfeedforward
import math,random
import xyplotter

def plote(plot,A,machine,weights):

   x1 = []
   y  = []
   y1 = []
   for i in range(-19,20):
      x = i/(19.0)
      e = A*x*x
      e1 = machine.evaluate([x],weights)[0]
      x1.append(x)
      y.append(e)
      y1.append(e1)
   plot.resetwindow(-1.0,0.0,1.0,0.3,"Spring Energies")
   plot.plot(x1,y,"blue")
   plot.plot(x1,y1,"red")

alpha = 0.005
A = 0.2

beta = 1.4
sigmoid      = lambda x: math.tanh(beta*x)
sigmoidp     = lambda x: beta*(1.0/math.cosh(beta*x))**2
sigmoidpp    = lambda x: -2.0*beta*beta*math.tanh(beta*x)*(1.0/math.sech(beta*x))**2
xmoid1   = lambda x: x
xmoidp1  = lambda x: 1.0
xmoidpp1 = lambda x: 0.0
xmoid2   = lambda x: x*x * 0.5
xmoidp2  = lambda x: x
xmoidpp2 = lambda x: 1.0
xmoid3   = lambda x: x*x*x * (1.0/6.0) 
xmoidp3  = lambda x: 3*x*x * (1.0/6.0)
xmoidpp3 = lambda x: 6*x   * (1.0/6.0)

xmoid4   = lambda x: x*x*x*x * (1.0/24.0)
xmoidp4  = lambda x: 4*x*x*x * (1.0/24.0)
xmoidpp4 = lambda x: 12*x*x *  (1/0/24.0)

#bias = [[0.01],[0.01],[0.001],[0.0001],[0.00001],[0.0000001]]
bias = [[0.01],[0.01],[0.001]]
bias = []
machine = myfeedforward.MyFeedForward([1,1,1,1,1,1,1,1],[sigmoid,sigmoid,sigmoid,xmoid3,sigmoid,xmoid2,sigmoid,xmoid1],[sigmoidp,sigmoidp,sigmoidp,xmoidp3,sigmoidp,xmoidp2,sigmoidp,xmoidp1],[sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp3,sigmoidpp,xmoidpp2,sigmoidpp,xmoidpp1],bias)

#machine = myfeedforward.MyFeedForward([1,1,1,1,1,1,1],[sigmoid,sigmoid,sigmoid,sigmoid,xmoid3,xmoid2,xmoid1],[sigmoidp,sigmoidp,sigmoidp,sigmoidp,xmoidp3,xmoidp2,xmoidp1],[sigmoidpp,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp3,xmoidpp2,xmoidpp1],bias)

#machine = myfeedforward.MyFeedForward([1,1,1,1,1,1],[sigmoid,sigmoid,sigmoid,sigmoid,xmoid2,xmoid1],[sigmoidp,sigmoidp,sigmoidp,sigmoidp,xmoidp2,xmoidp1],[sigmoidpp,sigmoidpp,sigmoidpp,sigmoidpp,xmoidpp2,xmoidpp1],bias)
weights = machine.initial_w()

for i in range(len(weights)):
   weights[i] = 1.0


plot  = xyplotter.xyplotter(-1.0,0.0,1.0,1.0, "Simple Plot for Raymond",2)

plote(plot,A,machine,weights)


for i in range(1000000):
   x = 2.5*random.random()-1.25
   e = A*x*x
   e1 = machine.evaluate([x],weights)[0]
   gg = machine.w_energy_gradient([x],[e],weights)
   error = gg[0]
   g1 = gg[1]
   for j in range(len(g1)):
      weights[j] -= alpha*g1[j]
   #print "x,e,e1,error,w=",x,e,e1,error,weights
   print x,e,e1,(e1-e)**2,error,weights[0],weights[1],weights[2],weights[3]

   if ((i%10000)==0):
      plote(plot,A,machine,weights)


