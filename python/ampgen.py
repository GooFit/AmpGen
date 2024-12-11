from ctypes import cdll, c_double, c_int, c_char_p, c_uint
from numpy import ctypeslib
import numpy as np
from tempfile import mkstemp
import subprocess 

def make_events(eventType, nEvts):
  lib = cdll.LoadLibrary("libAmpGen.so")
  gen_phsp = self.lib.python__Generate
  gen_phsp.argtypes = [ c_char_p, np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='C'), c_uint]
  nd = len(eventType.split(' ')) -1 ## number of final state particles
  c_evtType = c_char_p(eventType.encode('utf-8'))
  evts = np.zeros( shape = ( nEvts, 4*nd ) )
  self.gen_phsp( c_evtType, evts, nEvts )
  return evts 

def get_functions( libname ): 
    cmd = ["nm", "-D", libname ] 
    result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True) 
    lines = result.stdout.split('\n')
    all_functions = []
    for line in lines : 
        tokens = line.split(' ')
        if len(tokens) == 3 and tokens[1] == "T"  : 
            all_functions.append( tokens[2]  )

    return all_functions

class AmpGenModel:

    def __init__(self, libname):
        self.double_etype = c_double*16  ### FIXME 
        self.lib = cdll.LoadLibrary(libname)
        self.functionList = get_functions(libname)
        self.lib.FCN.argtypes = [self.double_etype, c_int] # E, parity
        self.lib.FCN.restype = c_double

    def FCN(self, event,*, cp=+1): 
        return self.lib.FCN(self.double_etype(*event), cp) 
        
def make_ampgen_model( options_file ) :
    fid,src=mkstemp( suffix='.cpp')
    print(src)
    lib=src.replace(".cpp",".so")
    output_log=src.replace(".cpp",".log")
    output_log_file=open(output_log, 'w')
    cmd = ["AmpGen", options_file, f"--Output {src}"]
    result = subprocess.run(cmd, stdout=output_log_file, stderr=output_log_file, check=False, text=True)
    return AmpGenModel( lib )
