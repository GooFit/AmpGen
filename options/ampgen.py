from ctypes import cdll, c_double, c_int, c_char_p, c_uint, POINTER, c_void_p
from numpy import ctypeslib
import numpy as np
from collections import namedtuple

# Compile with:
# clang++ -std=c++14 -shared -rdynamic -Ofast -fPIC model2.cpp -o model2.so

MatrixElement = namedtuple('MatrixElement',('pname', 'name', 'vars', 'amp', 'coef', 'coef_nopol'))

class Generator:
    __slots__ = ('lib','gen_phsp')

    def __init__(self):
        self.lib = cdll.LoadLibrary("libAmpGen.so")
        self.gen_phsp = self.lib.PyGenerate
        self.gen_phsp.argtypes = [ c_char_p
        , np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='C')
        , c_uint]
    
    def makeEvents(self, eventType, nd, nEvts):
        c_evtType = c_char_p(eventType.encode('utf-8'))
        evts = np.zeros( shape = ( nEvts, 4*nd ) )
        self.gen_phsp( c_evtType, evts, nEvts )
        return evts 

class FixedLib:
    __slots__ = ('double_etype', 'lib', 'matrix_elements', 'normalization')

    def __init__(self, libname):
        self.double_etype = c_double*16
        lib = cdll.LoadLibrary(libname)

        lib.matrix_elements_n.restype = c_int
        n_matel = lib.matrix_elements_n()

        lib.normalization.restype = c_double
        self.normalization = lib.normalization()

        lib.matrix_elements.argtypes = [c_int] # Number of mat ele
        lib.matrix_elements.restype = c_char_p
        progNames = [lib.matrix_elements(n).decode('ascii') for n in range(n_matel)]
        
        names_f = [getattr(lib, h+'_name') for h in progNames]
        for name_f in names_f: name_f.restype = c_char_p
        
        names = [name_f().decode('ascii') for name_f in names_f]

        plen_f = [getattr(lib, h+'_pSize') for h in progNames]
        for l in plen_f: l.restype = c_int
        plen = [l() for l in plen_f]

        pval_f = [getattr(lib, h+'_pVal') for h in progNames]
        for l in pval_f: l.restype = c_double
        pval = [[l(n) for n in range(q)] for q,l in zip(plen,pval_f)]

        pfuc_f = [getattr(lib, h ) for h in progNames]
        pfunc = []
        for l,f,h in zip(plen, pfuc_f, progNames):
            f.argtypes = [POINTER(c_double), POINTER(c_double), self.double_etype, c_double*l]
            def wrapper(P,E,*,l=l,f=f):
                real = c_double()
                imag = c_double()
                double_v = c_double * l
                f(real, imag, self.double_etype(*P), double_v(*E))
                return complex(real.value, imag.value)
            wrapper.__qualname__ = h
            wrapper.__name__     = h 
            pfunc.append(wrapper)

        lib.coefficients.argtypes = [c_int]*3 # n, which, parity
        lib.coefficients.restype = c_double
        parity   = [complex(lib.coefficients(n,0,-1),lib.coefficients(n,1,-1)) for n in range(n_matel)]
        noparity = [complex(lib.coefficients(n,0,+1),lib.coefficients(n,1,+1)) for n in range(n_matel)]
        self.matrix_elements = [MatrixElement(*l) for l in  zip(progNames, names, pval, pfunc, parity, noparity)]
        self.lib = lib

    def __len__(self):
        return len(self.matrix_elements)

    def __getitem__(self, item):
        return self.matrix_elements[item]

    def collect_vars(self):
        return {n:v for me in self.matrix_elements for n,v in zip(me.varnames, me.vars)}

    def FCN(self, PS, *, parity=1):
        self.lib.FCN.argtypes = [self.double_etype, c_int] # E, parity
        self.lib.FCN.restype = c_double
        return self.lib.FCN(self.double_etype(*PS), parity)

    def FCN_all(self, matrix, *, amps=None, parity=1, mt=True):
        FCN = self.lib.FCN_mt if mt else self.lib.FCN_all
        matrix = np.asarray(matrix, dtype=np.double, order='C')
        events = matrix.shape[0]
        FCN.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='C'),
            np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='C'),
            c_uint, c_int, POINTER(c_double) ]
        out = np.empty((events,), dtype=np.double)
        if amps is not None:
            amps = ctypeslib.as_ctypes(np.as_array(amps, dtype=np.double))
            assert amps.shape == (len(self.matrix_elements)*2,)
        print( out.shape )
        FCN(out, matrix, events, parity, amps)
        return out
