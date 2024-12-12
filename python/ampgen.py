from ctypes import cdll, c_double, c_int, c_char_p, c_uint, pointer, addressof 
from numpy import ctypeslib
import numpy as np
from tempfile import mkstemp
import subprocess 
import os 

def make_events(eventType, nEvts):
  lib = cdll.LoadLibrary("libAmpGen.so")
  gen_phsp = lib.python__generate
  gen_phsp.argtypes = [ c_char_p, np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='C'), c_uint]
  nd = len(eventType.split(' ')) -1 ## number of final state particles
  c_evtType = c_char_p(eventType.encode('utf-8'))
  evts = np.zeros( shape = ( nEvts, 4*nd ) )
  gen_phsp( c_evtType, evts, nEvts )
  return evts 

def get_matrix_element_dim(eventType): 
    lib = cdll.LoadLibrary("libAmpGen.so")
    dim_eventType = lib.python__EventType__dim
    dim_eventType.restype = c_int 
    rt = dim_eventType( c_char_p(eventType.encode('utf-8')) )
    print( rt ) 
    return ( 0xFFFF & ( rt >> 16 ) , 0xFFFF & rt  )


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
        self.lib = cdll.LoadLibrary(libname)
        self.lib.EventType.restype = c_char_p 
        particles = self.lib.EventType().decode('utf-8')
        self.mother        = particles.split(' ')[0]
        self.decayProducts = particles.split(' ')[1:]
        self.nParticles    = len(self.decayProducts) 
        self.functionList       = get_functions(libname)
        self.double_etype       = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
        self.lib.FCN.argtypes   = [self.double_etype, c_int] 

        self.lib.FCN.restype    = c_double
        self.matrix_element_dim = get_matrix_element_dim( particles ) 
        self.lib.all_amplitudes__py.argtypes = [self.double_etype, np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS') ] 
            
        self.lib.MatrixElements.restype = c_char_p 
        self.labels = self.lib.MatrixElements().decode('utf-8').split('\n')
        self.nAmplitudes = len(self.labels)
    def FCN(self, event, cp=+1):  return self.lib.FCN(event, cp) 
    
    def eventTypeString(self):
        return self.lib.EventType().decode("utf-8")
    
    def phsp_event(self, batch=16, which=0): 
        return make_events(self.eventTypeString(), batch)[which]

    def print_amplitude(self, event): 
        dim = self.matrix_element_dim[0] * self.matrix_element_dim[1]
        amp = np.zeros( shape = (dim * self.nAmplitudes, 2) )
        self.lib.all_amplitudes__py(event , amp)
        for me in range( 0, len(self.labels)): 
            elements = ""
            for i in range(0, dim) : 
                elements = f"{elements} ({amp[dim*me+i][0]}, {amp[dim*me+i][1]}) "
            print( f"{self.labels[me]} {elements}" )
        
    def amplitudes(self, event): 
        dim = self.matrix_element_dim[0] * self.matrix_element_dim[1]
        amp = np.zeros( shape = (self.matrix_element_dim[0] * self.matrix_element_dim[1] * len(self.labels), 2) )
        self.lib.all_amplitudes__py(event , amp)
        rt = dict() 
        for me in range( 0, len(self.labels)): 
            elements = []
            for i in range(0, dim) : 
                elements.append( amp[dim*me+i][0] ) # real part 
                elements.append( amp[dim*me+i][1] ) # imaginary part 
            rt[ self.labels[me] ] = elements 
        return rt  

def make_ampgen_model( options_file, src = "" ) :
    if src == "" : src=mkstemp( suffix='.cpp')[1]
    lib=src.replace(".cpp",".so")
    if not os.path.isfile(lib) :
        cmd = ["AmpGen", options_file, f"--Output {src}", "--Normalise=0", "--CompilerWrapper::Verbose=1"]
        print( f"Running command {cmd}")
        subprocess.run(cmd, check=True, text=True)
    return AmpGenModel(os.path.abspath(lib))

