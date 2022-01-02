import turtle


def _line(t,x1,y1,x2,y2):
   t.penup()
   t.goto(x1,y1)
   t.pendown()
   t.goto(x2,y2)

def _fillsquare(t,x1,y1,x2,y2):
   t.penup()
   t.goto(x1,y1)
   t.color("black","white")
   t.pendown()
   #t.fill(True)
   t.begin_fill()
   t.goto(x1,y2)
   t.goto(x2,y2)
   t.goto(x2,y1)
   t.goto(x1,y1)
   #t.fill(False)
   t.end_fill()

   t.penup()
   t.goto(x1, y1+0.25*(y2-y1))
   t.pendown()
   t.goto(x2, y1+0.25*(y2-y1))

   t.penup()
   t.goto(x1,y1+0.5*(y2-y1))
   t.pendown()
   t.goto(x2,y1+0.5*(y2-y1))

   t.penup()
   t.goto(x1, y1+0.75*(y2-y1))
   t.pendown()
   t.goto(x2, y1+0.75*(y2-y1))

   t.penup()
   t.goto(x1 + 0.25*(x2-x1),y1)
   t.pendown()
   t.goto(x1 + 0.25*(x2-x1),y2)

   t.penup()
   t.goto(x1 + 0.5*(x2-x1),y1)
   t.pendown()
   t.goto(x1 + 0.5*(x2-x1),y2)

   t.penup()
   t.goto(x1 + 0.75*(x2-x1),y1)
   t.pendown()
   t.goto(x1 + 0.75*(x2-x1),y2)
   t.penup()

def _labels(t,x1,y1,x2,y2,dx,dy,ltype=None):
   if (ltype==None): ltype=1
   if (ltype==2):
      xname = "Path"
      yname = "Relative Energy"
      yname2 = "kcal/mol"
      xf  = "%.2f"
      yf  = "%.2f"
      dxf = 0.40
      fnt = ('Arial', 12, 'normal')
      t.penup()
      t.goto(0.5*(x1+x2)-0.1*dx,y1-0.90*dy)
      t.pendown()
      t.write("%s" % (xname), font=fnt)
      for i in range(len(yname)):
         t.penup()
         t.goto(x1-0.80*dx,y2-(2.70+i*0.31)*dy)
         t.pendown()
         t.write("%s" % (yname[i]), font=fnt)
      for i in range(len(yname2)):
         t.penup()
         t.goto(x1-0.6*dx,y2-(3.85+i*0.31)*dy)
         t.pendown()
         t.write("%s" % (yname2[i]), font=fnt)
   elif (ltype==3):
      xname = "Path"
      yname = "Residual Error a.u."
      xf  = "%.2f"
      yf  = "%.6f"
      dxf = 0.60
      fnt = ('Arial', 12, 'normal')
      t.penup()
      t.goto(0.5*(x1+x2)-0.1*dx,y1-0.90*dy)
      t.pendown()
      t.write("%s" % (xname), font=fnt)
      for i in range(len(yname)):
         t.penup()
         t.goto(x1-0.80*dx,y2-(2.70+i*0.31)*dy)
         t.pendown()
         t.write("%s" % (yname[i]), font=fnt)
   elif (ltype==4):
      xname = "Distance a.u."
      yname = "Energy a.u."
      xf  = "%.2f"
      yf  = "%.6f"
      dxf = 0.60
      fnt = ('Arial', 12, 'normal')
      t.penup()
      t.goto(0.5*(x1+x2)-0.1*dx,y1-0.90*dy)
      t.pendown()
      t.write("%s" % (xname), font=fnt)
      for i in range(len(yname)):
         t.penup()
         t.goto(x1-0.80*dx,y2-(2.70+i*0.31)*dy)
         t.pendown()
         t.write("%s" % (yname[i]), font=fnt)
   else:
      xf  = "%.1e"
      yf  = "%.9e"
      dxf = 0.98
      fnt = ('Arial', 8, 'normal')
     

   t.penup()
   t.goto(x1,y1-0.5*dy)
   t.pendown()
   t.write(xf % (x1), font=fnt)

   t.penup()
   t.goto(x1+0.25*(x2-x1),y1-0.5*dy)
   t.pendown()
   t.write(xf % (x1+0.25*(x2-x1)), font=fnt)

   t.penup()
   t.goto(0.5*(x1+x2),y1-0.5*dy)
   t.pendown()
   t.write(xf % (0.5*(x1+x2)), font=fnt)

   t.penup()
   t.goto(x1+0.75*(x2-x1),y1-0.5*dy)
   t.pendown()
   t.write(xf % (x1+0.75*(x2-x1)), font=fnt)

   t.penup()
   t.goto(x2,y1-0.5*dy)
   t.pendown()
   t.write(xf % (x2), font=fnt)

   t.penup()
   t.goto(x1-dxf*dx,y1-0.12*dy)
   t.pendown()
   t.write(yf % (y1), font=fnt)

   t.penup()
   t.goto(x1-dxf*dx,y1+0.25*(y2-y1)-0.12*dy)
   t.pendown()
   t.write(yf % (y1+0.25*(y2-y1)), font=fnt)

   t.penup()
   t.goto(x1-dxf*dx,0.5*(y1+y2)-0.12*dy)
   t.pendown()
   t.write(yf % (0.5*(y1+y2)), font=fnt)

   t.penup()
   t.goto(x1-dxf*dx,y1+0.75*(y2-y1)-0.12*dy)
   t.pendown()
   t.write(yf % (y1+0.75*(y2-y1)), font=fnt)

   t.penup()
   t.goto(x1-dxf*dx,y2-0.12*dy)
   t.pendown()
   t.write(yf % (y2), font=fnt)
   t.penup()



class xyplotter:
   def __init__(self,xmin,ymin,xmax,ymax,title=None,ltype=None):
      dx = 0.25*(xmax-xmin)
      dy = 0.10*(ymax-ymin)
      self.xmin = xmin
      self.xmax = xmax
      self.ymin = ymin
      self.ymax = ymax
      self.dx   = dx
      self.dy   = dy
      self.title1 = title
      self.ltype  = ltype
      self.root = turtle.TK.Tk()
      if (title==None):
         self.root.title("xyplotter")
      else:
         self.root.title(title)
      #self.cv1  = turtle.TK.Canvas(self.root,width=750,height=600,bg="#ddffff")
      self.cv1  = turtle.TK.Canvas(self.root,width=500,height=400,bg="#ddffff")
      self.cv1.pack()
      self.s1 = turtle.TurtleScreen(self.cv1)
      self.s1.setworldcoordinates(xmin-dx,ymin-dy,xmax+dx,ymax+dy)
      self.t1 = turtle.RawTurtle(self.s1)
      self.t1.speed("fastest")
      self.t1.screen.tracer(10000)
      self.t1.hideturtle()
      self.t1.pensize(1)
      self.t1.color("black")
      _fillsquare(self.t1,self.xmin,self.ymin,self.xmax,self.ymax)
      _labels(self.t1,self.xmin,self.ymin,self.xmax,self.ymax,self.dx,self.dy,self.ltype)
      self.s1.update()

   def get_xmax(self): return self.xmax
   def get_xmin(self): return self.xmin
   def get_ymax(self): return self.ymax
   def get_ymin(self): return self.ymin
   def resetwindow(self,xmin,ymin,xmax,ymax,title=None):
      dx = 0.25*(xmax-xmin)
      dy = 0.10*(ymax-ymin)
      self.xmin = xmin
      self.xmax = xmax
      self.ymin = ymin
      self.ymax = ymax
      self.dx   = dx
      self.dy   = dy
      self.title1 = title
      self.s1 = turtle.TurtleScreen(self.cv1)
      self.s1.setworldcoordinates(xmin-dx,ymin-dy,xmax+dx,ymax+dy)
      self.t1 = turtle.RawTurtle(self.s1)
      self.t1.speed("fastest")
      self.t1.screen.tracer(10000)
      self.t1.hideturtle()
      self.t1.pensize(1)
      self.t1.color("black")
      _fillsquare(self.t1,self.xmin,self.ymin,self.xmax,self.ymax)
      _labels(self.t1,self.xmin,self.ymin,self.xmax,self.ymax,self.dx,self.dy,self.ltype)
      self.s1.update()
   
   def title(self,msg):
      self.t1.penup()
      self.t1.goto(0.5*(self.xmin+self.xmax)-0.08*(len(msg))*self.dx,self.ymax+0.4*(self.dy))
      self.t1.pendown()
      self.t1.write("%s" % msg,font=('Arial', 12, 'normal'))
      self.t1.penup()
      self.s1.update()

   def reset(self):
      self.t1.screen.tracer(0)
      self.t1.hideturtle()
      _fillsquare(self.t1,self.xmin,self.ymin,self.xmax,self.ymax)
      self.s1.update()

   def plot(self,x,y,color=None):
      n = len(x)
      if (n>1):
         if (color==None):
            self.t1.color("black")
         else:
            self.t1.color(color)
         self.t1.screen.tracer(10000)
         self.t1.hideturtle()
         self.t1.penup()
         self.t1.goto(x[0],y[0])
         self.t1.pendown()
         for i in range(1,n):
            self.t1.goto(x[i],y[i])
         self.s1.update()

   def plot1(self,y,color=None):
      n = len(y)
      if (n>1):
         if (color==None):
            self.t1.color("black")
         else:
            self.t1.color(color)
         self.t1.screen.tracer(10000)
         self.t1.hideturtle()
         self.t1.penup()
         x = 0.0
         self.t1.goto(x,y[0])
         self.t1.pendown()
         for i in range(1,n):
            x = float(i)/float(n-1)
            self.t1.goto(x,y[i])
         self.s1.update()

   def dotplot1(self,y,color=None):
      n = len(y)
      if (n>0):
         if (color==None):
            ctmp = "black"
         else:
            ctmp = color
         self.t1.screen.tracer(10000)
         self.t1.hideturtle()
         self.t1.penup()
         for i in range(n):
            x = float(i)/float(n-1)
            self.t1.goto(x,y[i])
            self.t1.dot(5,ctmp)
         self.s1.update()


   def dotplot(self,x,y,color=None):
      n = len(x)
      if (n>0):
         if (color==None):
            ctmp = "black"
         else:
            ctmp = color
         self.t1.screen.tracer(10000)
         self.t1.hideturtle()
         self.t1.penup()
         for i in range(0,n):
            self.t1.goto(x[i],y[i])
            self.t1.dot(5,ctmp)
         self.s1.update()

   def dotplot2(self,xy1,color1=None,xy2=None,color2=None):
      n1 = len(xy1)
      if (color1==None): ctmp1 = "black"
      else:              ctmp1 = color1
      if (xy2!=None): n2 = len(xy2)
      else:           n2 = 0
      if (color2==None): ctmp2 = "black"
      else:              ctmp2 = color2

      self.t1.speed("fastest")
      self.t1.screen.tracer(0)
      self.t1.hideturtle()
      _fillsquare(self.t1,self.xmin,self.ymin,self.xmax,self.ymax)
      if (n1>0):
         for i in range(n1):
            self.t1.goto(xy1[i][0],xy1[i][1])
            self.t1.dot(2,ctmp1)
      if (n2>0):
         for i in range(n2):
            self.t1.goto(xy2[i][0],xy2[i][1])
            self.t1.dot(2,ctmp2)
      self.s1.update()

   def dotplot2b(self,xy1,color1=None,xy2=None,color2=None):
      n1 = len(xy1)
      if (color1==None): ctmp1 = "black"
      else:              ctmp1 = color1
      if (xy2!=None): n2 = len(xy2)
      else:           n2 = 0
      if (color2==None): ctmp2 = "black"
      else:              ctmp2 = color2

      #self.t1.speed("fastest")
      self.t1.screen.tracer(10000)
      #self.t1.hideturtle()
      
      if (n1>0):
         for i in range(n1):
            self.t1.goto(xy1[i][0],xy1[i][1])
            self.t1.dot(1,ctmp1)
      if (n2>0):
         for i in range(n2):
            self.t1.goto(xy2[i][0],xy2[i][1])
            self.t1.dot(1,ctmp2)

   def update(self):
      self.s1.update()

   def print1(self,filename):
      self.cv1.postscript(file = filename, colormode='color')



