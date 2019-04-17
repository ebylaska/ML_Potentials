



#0.5*K*(0.125)**2 = 0.003
#dE = 1.170 - 1.167 = 0.003

K   = 2.0*0.003/(0.125)**2
E0  = -1.1705
r0 = 1.25 + 0.25/2.0
for i in range(20000):
   ii = i%500
   r = 1.25 + 0.25*ii/499
   E = E0 + 0.5*K*(r-r0)**2
   print "2"
   print E
   print "20.00000   0.00000   0.00000"
   print " 0.00000  20.00000   0.00000"
   print " 0.00000   0.00000  20.00000"
   print "1 H H %f %f %f %f %f %f" % (0.0,0.0,-0.5*r, 0.0,0.0,-K*(r-r0))
   print "2 H H %f %f %f %f %f %f" % (0.0,0.0, 0.5*r, 0.0,0.0, K*(r-r0))

