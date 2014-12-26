print(''' Module for manipulating nus schedule lists in the bruker and omega formats ''')
print(''' To use: nus.Omega2Bruk('nuslist1','nuslist2) will output a bruker format nus list 
for use with hmsIST''')

import numpy as np
import pylab as pl

def Bruk2Omega1d(nuslist):
  x = open('nuslist','r')
  lines = x.readlines()
  o = open('increment_hsqc','w')
  o.write('increment_hsqc[0]={\n')
  length = len(lines)-1
  for n,line in enumerate(lines):
     line=line.strip('\n')
     if n!=length:
       print('%s,'%line)
       o.write('%s,\n'%line)
     else:
       print('%s'%line)
       o.write('%s\n'%line)
  o.write('};')

def Omega2Bruk1d(nuslist):
  x = open(nuslist,'r')
  linesx = x.readlines()
  x.close()
  linesx=linesx[1:-1]
  
  array=[]
  for line in linesx:
    cols = line.split(',')
    try:
      for i in cols:
        array.append(int(i))
    except ValueError:
      pass
  print(array)
  o = open('nuslist','w')

  for i in array:
    o.write('%d\n'%(i))


def Omega2Bruk(nuslist1,nuslist2):
  x = open(nuslist1,'r')
  linesx = x.readlines()
  x.close()
  linesx=linesx[1:-1]
  y = open(nuslist2,'r')
  linesy = y.readlines()
  y.close()
  linesy=linesy[1:-1]

  array1=[]
  for line in linesx:
    cols = line.split(',')
    try:
      for i in cols:
        array1.append(int(i))
    except ValueError:
      pass
  print(array1)

  array2=[]
  for line in linesy:
    cols = line.split(',')
    try:
      for i in cols:
        array2.append(int(i))
    except ValueError:
      pass
  print(array2)

  o1 = open('sched2d','w')

  for i,j in zip(array1,array2):
    o1.write('%d %d\n'%(i,j))


def Bruk2Omega(nuslist):
  x = open(nuslist,'r')
  lines = x.readlines()
  x.close()

  file = str(nuslist)

  o1 = open('increment1','w')
  o2 = open('increment2','w')

  #o1.write(file+'256 [0] = {\n')
  o1.write('increment1 [0] = {\n')
  o2.write('increment2 [0] = {\n')
  ''' Calculating length of list so the last point is formatted correctly '''
  length = len(lines)-1
  for n,line in enumerate(lines):
    cols = line.split()
    inc1 = int(cols[0])
    inc2 = int(cols[1])
    if n!=length:
      o1.write('%d,\n'%inc1)
      o2.write('%d,\n'%inc2)
    else: 
      o1.write('%d\n'%inc1)
      o2.write('%d\n'%inc2)
      
  o1.write('};')
  o2.write('};')

def PlotSched(nuslist):
  a = np.genfromtxt(nuslist)
  t1= a[:,0]
  t2= a[:,1]
  print(len(t1))
  print(len(t2))
  
  fig=pl.figure()
  ax =fig.add_subplot(111)
  ax.scatter(t1,t2)
  ax.set_xlabel('t1')
  ax.set_ylabel('t2')
  print('saved plot as test.pdf')
  pl.savefig('test.pdf')  

def MakeLatexTable(nuslist):
  sched = open(nuslist,'r')
  x = sched.readlines()
  sched.close()
  for line in x:
    cols = line.split()
    print('%d & %d \\\ '%(int(cols[0]),int(cols[1])))
#def PlotSchedule(nuslist):
