import numpy as np

def sliding_window(sequence,window=19):
  seqs = [sequence[i:i+window] for i in range(len(sequence)-(window-1))]
  res  = np.arange((window/2+1),len(seqs)+(window/2+1))
  #print(res,seqs)
  return res,seqs

def sliding_window_function(sequence,function,window=19,**kwargs):
  ''' Performs function over each element of the sliding window '''
  res,seqs = sliding_window(sequence,window)
  if kwargs:
    result_dict  = dict([(n+1,function(i=i,**kwargs)) for n,i in enumerate(seqs)])
  else:
    result_dict  = dict([(n+1,function(i=i)) for n,i in enumerate(seqs)])
  return result_dict

def plot_over_sliding_window(i,x,y,z):
  return i+x+y+z
i=1
def f3(i):
  pass
def f1(i,**args):
  return args
def f2(i,**kwargs):
  return f1(i,**kwargs)




if __name__ == "__main__":
  seq = "ASDFGHGSHGSAHDGAHSDGASHDGAHSDGAHSDGHASGDASD"
  res,seqs = sliding_window(seq)
  
#  result = sliding_window_function(sequence=seq,window=19,function=plot_over_sliding_window,x="A",y="B",z="C")
  result = sliding_window_function(sequencs=seq,window=19,function=f3)
  print seqs
  print result

#print f1(i,a="a",b="b")
#print f2(i,a="a",b="b")
