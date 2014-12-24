import matplotlib.pyplot as pl
import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter
from collections import defaultdict
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

''' My modules  '''
from .BioTools import hmoment,plotwheel

def Data_URI_mpl(format="png"):
  import io
  import base64
  ''' Converts matplot lib image into a base64 string that can be
      rendered in a browser using the following code:
      <img src="data:image/png;base64,{{ plot }}" class="img-responsive" alt="name" style="width:auto;height:auto"> 
  ''' 
  sio = io.BytesIO()
  pl.savefig(sio, format=format)
  imgStr = base64.encodebytes(sio.getvalue()).decode()

  return imgStr


def plotWheel(seq,m,scale):
  r,theta,Hm,H_av,seq,direction = hmoment(clean(seq),m=m,scale=scale)

  f,ax = plot_wheel(r,theta,Hm,H_av,seq,direction,type=m)
  pl.annotate((r'H$\Phi$ = %8.3f   $\mu_H$ = %8.3f'%(H_av,Hm))
             ,(-np.pi/2,5.0),ha='center',va='center',fontsize=16)

  ''' writing data to base64 string '''
  return Data_URI_mpl(format="png")

def plotHydropathy(sequence,window,scale_choices):
  f  = pl.figure(facecolor='white')
  ax = f.add_subplot(211)
  ax2 = f.add_subplot(212)
  for i in scale_choices:
    scale = get_scale(i)
    data = hydrophobicMoment(sequence=str(sequence),hscale=scale,window=int(window))
    res,data = data.getHmoments()
    h = data[:,2]
    hm= data[:,3]
    ax.plot(res,h,label='%s, window=%d'%(i,window))
    ax2.plot(res,hm)
    axes = (ax,ax2)
    for axis in axes:
      axis.xaxis.set_minor_locator(minorLocator)
      axis.xaxis.set_major_locator(majorLocator)
    ax2.set_xlabel('$Residue$',fontsize=20)
    ax2.set_ylabel('$Hydrophobic$ $moment$',fontsize=20)
    ax.set_ylabel('$Hydrophobicity$',fontsize=20)
    leg=ax.legend(loc=1,bbox_to_anchor=(1.05,1.2),prop={'size':10},ncol=3)
    leg.draw_frame(False)

  ''' writing data to base64 string '''
  return Data_URI_mpl(format="png")
