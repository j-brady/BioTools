#!/usr/bin/env python3
''' Hydropathy Data '''

''' White-Whimley hydropathy scale with Interface,Octanol and Octanol-Interface values in deltaG(kcal/mol)  '''
WhiteWhimley= {'K': (0.99, 2.8, 1.81), 'V': (0.07, -0.46, -0.53),
               'F': (-1.13, -1.71, -0.58), 'C': (-0.24, -0.02, 0.22),
               'D': (1.23, 3.64, 2.41), 'E0': (-0.01, 0.11, 0.12),
               'M': (-0.23, -0.67, -0.44), 'L': (-0.56, -1.25, -0.69),
               'N': (0.42, 0.85, 0.43), 'Y': (-0.94, -0.71, 0.23),
               'I': (-0.31, -1.12, -0.81), 'Q': (0.58, 0.77, 0.19),
               'T': (0.14, 0.25, 0.11), 'H+': (0.96, 2.33, 1.37),
               'G': (0.01, 1.15, 1.14), 'D0': (-0.07, 0.43, 0.5),
               'W': (-1.85, -2.09, -0.24), 'E': (2.02, 3.63, 1.61),
               'S': (0.13, 0.46, 0.33), 'P': (0.45, 0.14, -0.31),
               'A': (0.17, 0.5, 0.33), 'R': (0.81, 1.81, 1.0),
               'H': (0.17, 0.11, -0.06)}

''' Eisenberg consensus scale '''
Eisenberg  =  {'I':0.73,'F':0.61,'V':0.54,'L':0.53,'W':0.37,
               'M':0.26,'A':0.25,'G':0.16,'C':0.04,'Y':0.02,
               'P':-0.07,'T':-0.18,'S':-0.26,'H':-0.40,'E':-0.62,
               'N':-0.64,'Q':-0.69,'D':-0.72,'K':-1.10,'R':-1.76}

''' FaucherePliska scale '''
FaucherePliska={'A':0.310,'R':-1.010,'N':-0.600,'D':-0.770,
                'C':1.540,'Q':-0.220,'E':-0.640,'G':0.000,
                'H':0.130,'I':1.800,'L':1.700,'K':-0.990,
                'M':1.230,'F':1.790,'P':0.720,'S':-0.040,  
                'T':0.260,'W':2.250,'Y':0.960,'V':1.220}

''' Kyte Doolittle scale'''    
KyteDoolittle={'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,
               'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
               'H':-3.2,'I':4.5,'L': 3.8,'K':-3.9,
               'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
               'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}

''' Hopp Woods '''
HoppWoods={'A':-0.5,'R':3.0,'N':0.2,'D':3.0,
                'C':2.5,'Q':0.2,'E':3.0,'G':0.0,
                'H':-0.5,'I':-1.8,'L':-1.8,'K':3.0,
                'M':-1.3,'F':-2.5,'P':0.0,'S':0.3,
                'T':-0.4,'W':-3.4,'Y':-2.3,'V':-1.5}

''' Scale from EMBOSS program '''
H_scale={'A':0.62,'B':-0.84,'C':0.29,'D':-0.90,'E':-0.74,'F':1.19,'G':0.48,'H':-0.40,
         'I':1.38,'J':1.18,'K':-1.50,'L':1.06,'M':0.64,'N':-0.78,'O':-0.68,'P':0.12,
         'Q':-0.85,'R':-2.53,'S':-0.18,'T':-0.05,'U':-0.68,'V':1.08,'W':0.81,'X':-0.68,
         'Y':0.26,'Z':-0.78}	

''' Dictionary of scales '''
scale_dict = {'KyteDoolittle':KyteDoolittle,
              'HoppWoods':HoppWoods,
              'Eisenberg':Eisenberg,
              'FaucherePliska':FaucherePliska,
              'WhiteWhimley':WhiteWhimley,
              'Emboss':H_scale}

def get_scale(scale):
  return scale_dict[scale]

''' Turn propensity  Monne et al. (1999) '''
propensity = {'I':0.5,'F':0.3,'V':0.4,'L':0.3,'W':0.9,
           'M':0.5,'A':0.4,'G':1.1,'C':0.7,'Y':0.9,
           'P':1.6,'T':0.7,'S':0.9,'H':1.4,'E':1.3,
           'N':1.6,'Q':1.4,'D':1.5,'K':1.4,'R':1.5}

normalised = {'I':0.9,'F':0.6,'V':0.7,'L':0.5,'W':1.6,
           'M':0.8,'A':0.7,'G':1.9,'C':1.1,'Y':1.5,
           'P':2.7,'T':1.2,'S':1.6,'H':2.3,'E':2.1,
           'N':2.6,'Q':2.3,'D':2.5,'K':2.3,'R':2.6}
