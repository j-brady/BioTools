#!/usr/bin/env python3
#import pylab as pl
import os
import glob
from collections import defaultdict

import matplotlib.pyplot as pl
import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np
from Bio import SeqIO

''' My modules '''
from .data_tables import scale_dict, get_scale


''' General functions '''
def clean(sequence):
    ''' Function cleans protein sequence
        removing spaces, tabs and carriage
        returns/newlines '''
    remove = ['\t','\n','\r',' ']
    for i in remove:
        sequence = sequence.replace(i,'')
    return sequence

def sliding_window(sequence,window=19):
  seqs = [sequence[i:i+window] for i in range(len(sequence)-(window-1))]
  res  = np.arange((window/2+1),len(seqs)+(window/2+1))
  #print(res,seqs)
  return res,seqs

def sliding_window_function(sequence,function,window=19):
  ''' Performs function over each element of the sliding window '''
  res,seqs = sliding_window(sequence,window)
  result_dict  = dict([(n+1,function(i)) for n,i in enumerate(seqs)])
  return result_dict

''' Calculate hmoment over sliding window '''
def sliding_window_hmoment(sequence,window):
  ''' Returns dict containing {resnu:[Hm,H_av,seq]} starting with resnu=1 for the first window'''
  res, seqs = sliding_window(sequence,window)
  HmH  = dict([(n+1,hmoment(i)[-4:-1]) for n,i in enumerate(seqs)])
  
  return HmH

def sum_scale(sequences,hscale):
  data = []
  for seq in sequences: 
    sum = np.sum(np.vstack([hscale[i] for i in seq]),axis=0)
    data.append(sum)
  return data

def average_scale(sequences,hscale):
  data = []
  for seq in sequences:
    mean = np.mean(np.vstack([hscale[i] for i in seq]),axis=0) 
    data.append(mean)
  return data

def getHmoment(sequence,hscale,d,m):
    r = np.array([1.5 for i in sequence])
    #angle = np.array([d*i for i,j in enumerate(seq)])
    #theta = angle/(360.0/(2 * np.pi))%(2 * np.pi)
    delta = (np.pi*2)/m
    seq = [i for i in sequence]
    theta =[delta*n for n,i in enumerate(seq)]
    #theta = theta[::-1]
    theta = np.array(theta)*-1
    H = np.array([hscale[i] for i in seq])
    H_av= np.mean(H)
    Hm = np.sqrt((np.square(np.sum(H*np.sin(theta))) +\
               (np.square(np.sum(H*np.cos(theta))) )))
    direction = np.arctan2(np.sum(H*np.sin(theta)),np.sum(H*np.cos(theta)))
    ''' Making max hmoment face down  '''
    theta = theta-(np.pi/2 + direction)
    direction = np.arctan2(np.sum(H*np.sin(theta)),np.sum(H*np.cos(theta)))
    Hm= Hm/len(seq)  
    return r,theta,H_av,Hm,seq,direction

def read_fasta(fasta):
  fasta = SeqIO.read(open(fasta,'rU'),'fasta')
  sequence = fasta.seq
  return sequence

''' Hydropathy calculation '''

class hydropathy(object):
  def __init__(self,sequence,hscale,window):
    self.sequence = sequence
    self.hscale   = get_scale(hscale)
    self.window   = window

  def sum_hydrophobicities(self):
    ''' If you are using the white whimley scale you get three plots '''
    if self.hscale is WhiteWhimley:
      seqs = [self.sequence[i:i+self.window] for i in range(len(self.sequence)-(self.window-1))]
      res  = np.arange((self.window/2+1),len(seqs)+(self.window/2+1))  
      data = []
      for seq in seqs: 
        sum = np.sum(np.vstack([[self.hscale[i][0],self.hscale[i][1],self.hscale[i][2]] for i in seq]),axis=0)
        data.append(sum)
        data = np.vstack(data)
        #HmH  = dict([(n+1,hmoment(i)[-4:-1]) for n,i in enumerate(seqs)])
      ''' Interface '''
      i   = data[:,0]
      ''' Octanol '''
      o   = data[:,1]
      ''' Octanol-Interface '''
      io  = data[:,2]

      plots=(('Interface Scale',i),('Octanol Scale',o),('Octanol-Interface Scale',io))
      return res,plots
    else:
      '''Sum or mean can be either "average_scale" or "sum_scale" function'''
      res,seqs = window_function(sequence=self.sequence,window=self.window)
      data     = sum_scale(sequences=seqs,hscale=self.hscale)
      return res,data
  
  def av_hydrophobicities(self):
    if self.hscale is WhiteWhimley:
      seqs = [self.sequence[i:i+self.window] for i in range(len(self.sequence)-(self.window-1))]
      res  = np.arange((self.window/2+1),len(seqs)+(self.window/2+1))  
      data = []
      for seq in seqs: 
        mean = np.mean(np.vstack([[self.hscale[i][0],self.hscale[i][1],self.hscale[i][2]] for i in seq]),axis=0)
        data.append(mean)
      data = np.vstack(data)
        #HmH  = dict([(n+1,hmoment(i)[-4:-1]) for n,i in enumerate(seqs)])
      ''' Interface '''
      i   = data[:,0]
      ''' Octanol '''
      o   = data[:,1]
      ''' Octanol-Interface '''
      io  = data[:,2]

      plots=(('Interface Scale',i),('Octanol Scale',o),('Octanol-Interface Scale',io))
      return res,plots
    else:
      '''Sum or mean can be either "average_scale" or "sum_scale" function'''
      res,seqs = window_function(self.sequence,window=self.window)
      data     = average_scale(sequences=seqs,hscale=self.hscale)
      return res,data
  
  def raw_hydrophobicities(self):  
    '''raw data over amino acids '''
    data = [self.hscale[i] for i in self.sequence]
    res  = np.arange(1,len(data)+1) 
    return res,data
  #def plot_results(self):
      
  def __str__(self):
    return 'You\'re sequence is: \n %s'%self.sequence


class hydrophobicMoment(hydropathy):
  '''  Note currently uses the value of m for period '''
  def __init__(self,sequence,hscale,window,d=100,m=3.6):
    hydropathy.__init__(self,sequence,hscale,window)
    self.sequence = clean(sequence)
    self.hscale   = get_scale(hscale)
    self.window   = window
    self.d        = 100 #Degrees per amino acid
    self.m        = 3.6 #Amimo acids per turn
    
  def getHmoments(self):
    res,seqs = window_function(sequence=self.sequence,window=self.window)
    data  = np.vstack(np.array([getHmoment(sequence=i,hscale=self.hscale,d=self.d,m=self.m) for i in seqs],dtype='object'))
    return res,data
  
  def getHmomentMax(self,res,data):
    ind_max = np.argmax(data[:,3])
    res_max = res[ind_max]
    h_max   = data[:,2][ind_max]
    hm_max  = data[:,3][ind_max]
    seq_max = data[:,4][ind_max]
    return res_max,seq_max,h_max,hm_max

############################################################
''' Helical Wheels '''

''' Colors for helical wheel plots '''
aa_face_color ={'A': '0.7', 'C': 'y', 'E': 'r', 'D': 'r',
                'G': '0.7', 'F': 'y', 'I': 'y', 'H': 'c',
                'K': 'b', 'M': 'y', 'L': 'y', 'N': 'orange',
                'Q': 'orange', 'P': 'purple', 'S': 'g', 'R': 'b',
                'T': 'g', 'W': 'y', 'V': 'y', 'Y': 'y'}

aa_text_color ={'A': 'k', 'C': 'k', 'E': 'w', 'D': 'w',
                'G': 'k', 'F': 'k', 'I': 'k', 'H': 'k',
                'K': 'w', 'M': 'k', 'L': 'k', 'N': 'k',
                'Q': 'k', 'P': 'k', 'S': 'w', 'R': 'w',
                'T': 'w', 'W': 'k', 'V': 'k', 'Y': 'k'}

def pastel(colour, weight=2.0):
    """ Convert colour into a nice pastel shade"""
    rgb = np.asarray(colorConverter.to_rgb(colour))
    # scale colour
    maxc = max(rgb)
    if maxc < 1.0 and maxc > 0:
        # scale colour
        scale = 1.0 / maxc
        rgb = rgb * scale
    # now decrease saturation
    total = rgb.sum()
    slack = 0
    for x in rgb:
        slack += 1.0 - x

    # want to increase weight from total to weight
    # pick x s.t.  slack * x == weight - total
    # x = (weight - total) / slack
    x = (weight - total) / slack

    rgb = [c + (x * (1.0-c)) for c in rgb]

    return rgb

def hmoment(seq,m='alpha',scale='Eisenberg'):
  helix_dict = {'alpha':3.6,'three_ten':3.0,'three_eleven':(11.0/3.0),'pi_helix':4.1}
  
  ''' Making compatible with more than 18 aa '''
  r = np.array([2.6 for i in seq])
  inds = np.array(range(len(r)),dtype='int')

  if m == 'alpha':
    factor = inds//int(18) # // performs int arithmetic
  elif m == 'pi_helix':
    factor = inds//int(41) # // performs int arithmetic
    r = r
  elif m == 'three_ten':
    factor = inds//int(3) # // performs int arithmetic
    r = r-2.
  elif m == 'three_eleven':
    factor = inds//int(11) # // performs int arithmetic
    r = r-1
    
  m = helix_dict[m]
  scale = get_scale(scale)

  r = r + (factor*0.90)
  
  '''  '''
  #angle = np.array([d*i for i,j in enumerate(seq)])
  #theta = angle/(360.0/(2 * np.pi))%(2 * np.pi)
  delta = (np.pi*2)/m
  seq = [i for i in seq]
  theta =[delta*n for n,i in enumerate(seq)]
  ''' for more than 18 '''
  #theta = theta[::-1]
  theta = np.array(theta)*-1
  print('theta') 
  print(theta)
  H = np.array([scale[i] for i in seq])
  H_av= np.mean(H)
  Hm = np.sqrt((np.square(np.sum(H*np.sin(theta))) +\
               (np.square(np.sum(H*np.cos(theta))) )))
  direction = np.arctan2(np.sum(H*np.sin(theta)),np.sum(H*np.cos(theta)))
  ''' Making max hmoment face down  '''
  theta = theta-(np.pi/2 + direction)
  direction = np.arctan2(np.sum(H*np.sin(theta)),np.sum(H*np.cos(theta)))
  Hm= Hm/len(seq)  
  return r,theta,Hm,H_av,seq,direction


def plot_wheel(r,theta,Hm,H_av,seq,direction,type):
  ''' Type can be 'alpha,pi_helix,three_ten or three_eleven' '''
#def plot_wheel(hmoment):
  theta = theta
  fig= pl.figure(facecolor='white')
  ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=True)
  print('this is r')
  print(r)
  ''' Only plotting lines for up to 18  aa'''
  if len(theta>18) and type=='alpha':
    ax.plot(theta[:18], r[:18], color='k',linewidth=1.5)
  elif len(theta>41) and type=='pi_helix':
    ax.plot(theta[:41], r[:41], color='k',linewidth=1.5)
  elif len(theta>3) and type=='three_ten':
    ax.plot(theta[:3], r[:3], color='k',linewidth=1.5)
  elif len(theta>11) and type=='three_eleven':
    ax.plot(theta[:11], r[:11], color='k',linewidth=1.5)
  else:
    ax.plot(theta,r,color='k',linewidth=1.5)

  ax.grid(False)
  ax.set_yticklabels([])
  ax.set_xticklabels([])
  print(direction)
  ax.set_frame_on(False)  
  for num,vals in enumerate(zip(theta,r,seq)):
    i,j,k = vals[0],vals[1],vals[2]

    if aa_face_color[k].startswith('0') or aa_face_color[k] is 'r' or aa_face_color[k] is 'b' or aa_face_color[k] is 'orange' or aa_face_color[k] is 'g':
      ax.plot(i, j, marker='o',markersize=30,mfc=aa_face_color[k],mec=aa_face_color[k], linewidth=3)
    else:
      ax.plot(i, j, marker='o',markersize=30,mfc=pastel(aa_face_color[k]),mec=pastel(aa_face_color[k]), linewidth=3)

    if type == 'three_ten':
      ax.set_rmax(7.0)
    elif type == 'three_eleven':
      ax.set_rmax(5.0)
    elif type == 'pi_helix':
      ax.set_rmax(5.0)
    else:
      ax.set_rmax(5.0)
#    print k
    ax.annotate(k,(i,j),color=aa_text_color[k],weight='bold',fontsize=20,ha='center',va='center')
  ax.annotate('N',(theta[0],r[0]+0.4),color='r',fontsize=12,ha='center',va='center')
  ax.annotate('C',(theta[-1],r[-1]+0.4),color='b',fontsize=12,ha='center',va='center')
  ''' Direction arrow '''
  ax.annotate('', (-np.pi/2,0.75),
                (-np.pi/2, 0),
                #xycoords="figure fraction", textcoords="figure fraction",
                ha="center", va="center",
                arrowprops=dict(arrowstyle="simple",
                                fc="w", ec="k",
                                ),
                )
  #ax.arrow((np.pi/2),0,(np.pi/2),0.75,width=0.1,arrowstyle="bar")
  ax.annotate(r'$\mu_H$',(direction,0.8),color='r',
              fontsize=20,va='center',ha='center')
              #arrowprops=dict(arrowstyle="->",
              #              connectionstyle="arc3"))
  return fig,ax

def plot_2dhist(HmH):
  ax = pl.subplot(111)
  #ax.scatter(x,y,color='pink')
  ax.set_xlabel('Hydrophobicity',fontsize=14)
  ax.set_ylabel('Hydrophobic Moment',fontsize=14)
  ax.set_ylim(0,0.8)
  ax.set_xlim(-0.55,0.7)

  pl.hexbin(x,y,gridsize=7,cmap=pl.cm.Reds)
  pl.colorbar(shrink=0.4)
  pl.hexbin(tm_h,tm_hm,gridsize=5,cmap=pl.cm.Blues)
  pl.colorbar(shrink=0.4,orientation = 'horizontal')
  pl.plot(-0.053,0.358,marker='*',markersize=20,color='k')
  #pl.tight_layout(h_pad = 1.5,w_pad=1.5)
  pl.savefig('splatter.pdf')
def annotate_wheel(ax):
  ''' Takes axis and dict of sequences '''
  name = entry['name']
  start= entry['start']
  end  = entry['end']
  ax.annotate('%s \n %d to %d'%(name,start,end),(0,0),ha='center',va='center',color='k')
  #ax.annotate('%d to %d'%(start[0],start[0]+window),(0,0),ha='center',va='center',color='k',fontsize=12)


''' Plotting Hm and H vs sequence  '''
def plot_hmh_vs_seq(HmH):
  fig = pl.figure()
  axis = fig.add_subplot(111,polar=False)
  data = np.array([[res,hm[0],hm[1]] for res,hm in sorted(HmH.iteritems())])
  res = data[:,0]
  hm  = data[:,1]
  h   = data[:,2]
  axis.plot(res,hm,'k',label='Hydrophobic moment')
  axis.plot(res,h,'r',label='Hydrophobicity')
  pl.legend(loc=1)
  pl.show()
#def get_res_hm_h(HmH):
    

def get_hm_max(HmH):
  data = np.array([[res,hm[0],hm[1],hm[2]] for res,hm in sorted(HmH.iteritems())],dtype='object')
  res = data[:,0]
  hm  = data[:,1]
  h   = data[:,2]
  seqs= data[:,3]
  max_index = np.argmax(hm)
  res_max,h_max,hm_max,seq_max = res[max_index],h[max_index],hm[max_index],seqs[max_index] 
  return res_max,h_max,hm_max,seq_max

def get_hm_distributions(list_of_HmH,name,sd=True):
    majorLocator   = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(1)
    ''' Returns dictionary of the sums of all the Hm values across the sequences
        i.e. res:sum_of_hm'''
    #data = reduce(lambda x, y: dict((k, v + y[k]) for k, v in x.iteritems()),list_of_HmH)
    ''' Sums values for each key over all the dictionaries '''
    total_number_of_sequences=len(list_of_HmH)
    output=defaultdict(float)
    std_dict = defaultdict(list)
    numbers=defaultdict(float)
    for d in list_of_HmH:
        for k,v in d.iteritems():
           output[k]+=v[0]
           std_dict[k].append(v[0])
           numbers[k]+=1.0

    averages = np.array([(v/numbers[k]) for k,v in sorted(output.iteritems())])
    res_sums = np.array([(k,v) for k,v in sorted(output.iteritems())])
    num_of_each = np.array([v for k,v in sorted(numbers.iteritems())])
    ''' Standard deviation over each residue '''
    sd = np.array([np.std(v) for k,v in sorted(std_dict.iteritems())])
    
    res = res_sums[:,0]
    sums= res_sums[:,1]
    ''' Scaled average '''
    gridsize=int(max(res))
    fig = pl.subplot(111)
    
    #fig.errorbar(res,averages,yerr=sd,fmt='.',ecolor='0.5',alpha='0.8',markerfacecolor='none',visible=False,zorder=1)
    ax=pl.hexbin(res,averages, C=num_of_each, gridsize=gridsize, cmap=pl.cm.jet, bins=None)
    fig.set_xlabel('Linker length')
    fig.set_ylabel('Average Hydrophobic moment')
    fig.set_ylim(0,0.65)
    cb = pl.colorbar()
    cb.set_label('number of sequences')
    fig.xaxis.set_major_locator(majorLocator)
    fig.xaxis.set_major_formatter(majorFormatter)
    fig.xaxis.set_minor_locator(minorLocator)
    
    
    pl.savefig('2dhist_%s.pdf'%name)
    pl.clf()
    ''' Plotting data '''
    fig = pl.subplot(211)
    fig2= pl.subplot(212)

    
    axes = (fig,fig2)
    for axis in axes:
      axis.xaxis.set_major_locator(majorLocator)
      axis.xaxis.set_major_formatter(majorFormatter)
      axis.xaxis.set_minor_locator(minorLocator)
    
    fig.bar(res,sums,width=0.75,align='center')
    fig.set_ylabel('Sum of hydrophobic moments')
    fig.set_xlabel('Linker length')
    if sd is True:
      fig2.bar(res,averages,width=0.75,align='center',yerr=sd)
    else:
      fig2.bar(res,averages,width=0.75,align='center')
    fig2.set_ylabel('Average of hydrophobic moments')
    fig2.set_xlabel('Linker length')
    pl.tight_layout()
    pl.savefig('SumsAvs%s.pdf'%name)
    pl.clf()
    return res,sums,averages,sd
    #return output,numbers



def plot_hm_vs_h(res,h,hm):
  ax = pl.subplot(111)
  ax.scatter(h,hm)
  ax.set_ylim(0,0.8)
  ax.set_xlim(-0.55,0.7)
  ax.set_xlabel('$Hydrophobicity$',fontsize=20)
  ax.set_ylabel('$Hydrophobic$ $Moment$',fontsize=20)
  #pl.savefig(name+'.pdf')
