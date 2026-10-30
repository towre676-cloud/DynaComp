"""Replayed finite-horizon guarantees for a deterministic C binary64 engine.

Separate multiplication and addition, sequential left-to-right accumulation,
no FMA/reassociation, RN-even and gradual underflow are the execution contract.
Compilation is explicit on first preparation; no compiled binaries are shipped.
"""
import atexit
import copy
import os
import platform
import sys
import ctypes
from fractions import Fraction as Q
import hashlib
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

from .machine import _certify_execution, _up, U, ETA, LIMIT
from .binary64 import from_bits
from .rational import q
from .validated import replay_validated

CONTRACT = ('finite stored binary64 inputs; sequential left-to-right dot from +0; '
            'separate RN-even product and addition; gradual underflow; no FMA')
_LIBRARY = None
_TEMP = None


def sequential_budget(magnitude, count):
    """2k rounded operations: gamma_(2k) M + 2k eta/(1-2ku)."""
    if type(count) is not int or count < 0 or 2*count*U >= 1:
        raise ValueError('Invalid dot length')
    if not magnitude:
        return Q(0)
    denominator = 1-2*count*U
    error = (2*count*U*magnitude+2*count*ETA)/denominator
    if magnitude+error > LIMIT:
        raise ValueError('Cannot exclude intermediate sequential overflow')
    return _up(error)


def certify_sequential(E, G, H, D, initial_center, initial_radius, input_radii,
                       tolerance, input_centers=None):
    return _certify_execution(E, G, H, D, initial_center, initial_radius,
                             input_radii, tolerance, input_centers,
                             sequential_budget, CONTRACT,
                             'dynacomp.sequential-binary64/1')


def replay_sequential(packet):
    if packet.get('schema') != 'dynacomp.sequential-binary64/1':
        raise ValueError('Unknown sequential execution schema')
    fresh = certify_sequential(*(packet[k] for k in ('E','G','H','D')),
                               packet['initial_center'], packet['initial_radius'],
                               packet['input_radii'], packet['tolerance'], packet['input_centers'])
    if fresh != packet:
        raise ValueError('Sequential execution certificate does not reconstruct')
    return True


def compose_sequential(source, tolerance):
    if source.get('schema') != 'dynacomp.validated-box/1':
        raise ValueError('A finite validated-box source is required')
    replay_validated(source)
    runtime = certify_sequential(*(source[k] for k in ('E','G','H','D')),
                                 source['initial_center'], source['initial_radius'],
                                 source['input_radii'], tolerance, source['input_centers'])
    combined = [[str(q(a)+q(b)) for a,b in zip(model,execution)] for model,execution in
                zip(source['certificate']['output_upper'], runtime['output_error_upper'])]
    worst = max(q(v) for row in combined for v in row)
    return {'schema':'dynacomp.sequential-composition/1','source':source,'execution':runtime,
            'tolerance':str(q(tolerance)), 'combined_output_upper':combined,
            'combined_upper':str(worst),
            'status':'tolerance_met' if worst <= q(tolerance) and source['reduced_contraction_certified'] is not False else 'not_certified',
            'scope':'finite sampled rational source/model plus specified sequential C binary64 execution'}


def replay_sequential_composition(packet):
    if packet.get('schema') != 'dynacomp.sequential-composition/1' or packet != compose_sequential(packet['source'],packet['tolerance']):
        raise ValueError('Sequential composition does not reconstruct')
    return True


def _compiler():
    """Select a compiler executable, never a shell command or user flags."""
    names = ([os.environ['DYNACOMP_CC']] if os.environ.get('DYNACOMP_CC') else
             (['x86_64-w64-mingw32-gcc', 'gcc'] if sys.platform == 'win32'
              else ['cc', 'gcc', 'clang']))
    for name in names:
        compiler = shutil.which(name)
        if not compiler:
            continue
        if sys.platform == 'win32':
            if platform.machine().lower() not in ('amd64', 'x86_64') or ctypes.sizeof(ctypes.c_void_p) != 8:
                raise RuntimeError('Native Windows execution requires x86-64 Python')
            target = subprocess.run([compiler, '-dumpmachine'], check=True,
                                    capture_output=True, text=True).stdout.strip().lower()
            if not (target.startswith('x86_64-') and 'mingw' in target):
                continue
        return compiler
    raise RuntimeError('Install a C compiler (Windows: x86-64 MinGW-w64 GCC); '
                       'DYNACOMP_CC may name its executable')


def _compile_command(compiler, source, target):
    flags = ['-std=c11', '-O2', '-shared', '-fno-fast-math',
             '-ffp-contract=off', '-frounding-math']
    if sys.platform != 'win32':
        flags += ['-fPIC']
    else:
        # SSE2 directly evaluates binary64; never use x87 excess precision.
        flags += ['-msse2', '-mfpmath=sse', '-static-libgcc']
    return [compiler, *flags, str(source), '-lm', '-o', str(target)]


def _cleanup_library(lib, folder):
    # Windows prevents deletion of loaded DLLs. Run before TemporaryDirectory's
    # finalizer, at process exit only (live PreparedSequential objects share it).
    if sys.platform == 'win32':
        from _ctypes import FreeLibrary
        FreeLibrary(lib._handle)
    folder.cleanup()


def _library():
    global _LIBRARY, _TEMP
    if _LIBRARY is not None:
        return _LIBRARY
    compiler = _compiler()
    source = Path(__file__).with_name('sequential.c')
    folder = tempfile.TemporaryDirectory(prefix='dynacomp-sequential-')
    target = Path(folder.name)/('runtime.dll' if sys.platform == 'win32' else 'runtime.so')
    command = _compile_command(compiler, source, target)
    started = time.perf_counter()
    subprocess.run(command,check=True,capture_output=True,text=True)
    build_seconds = time.perf_counter()-started
    lib = ctypes.CDLL(str(target))
    pointer = ctypes.POINTER(ctypes.c_double)
    lib.run.argtypes = [ctypes.c_size_t]*6+[pointer]*7
    lib.run.restype = ctypes.c_int
    lib.compiler_command = command[:-1]+['<temporary runtime>']
    lib.source_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
    lib.build_seconds = build_seconds
    lib.compiler_version = subprocess.run([compiler,'--version'],check=True,capture_output=True,text=True).stdout.splitlines()[0]
    atexit.register(_cleanup_library, lib, folder)
    _TEMP = folder
    _LIBRARY = lib
    return lib


class PreparedSequential:
    """Checked batched execution; array axes match PreparedValidated."""
    def __init__(self,packet):
        import numpy as np
        from .prepared import PreparedSwitching
        self._packet = copy.deepcopy(packet)
        if packet.get('schema') == 'dynacomp.sequential-composition/1':
            replay_sequential_composition(self._packet)
            execution = self._packet['execution']
        else:
            replay_sequential(self._packet)
            execution = self._packet
        if self._packet['status'] != 'tolerance_met':
            raise ValueError('Sequential plan has not met tolerance')
        self._lib = _library()
        self._e,self._g,self._h,self._d = [np.ascontiguousarray([[from_bits(v) for v in row]
                                                           for row in execution['loaded_words'][k]],dtype=np.float64)
                                       for k in ('E','G','H','D')]
        self._r,self._n,self._m,self._p = len(self._g),self._e.shape[1],self._h.shape[1],len(self._d)
        for value in (self._e,self._g,self._h,self._d):
            value.flags.writeable = False
        self._horizon = len(execution['input_radii'])
        def endpoints(center,radius):
            return [PreparedSwitching._float_endpoints(q(c)-q(r),q(c)+q(r)) for c,r in zip(center,radius)]
        initial = endpoints(execution['initial_center'],execution['initial_radius'])
        self._lower = np.asarray([lo for lo,hi in initial])[:,None]
        self._upper = np.asarray([hi for lo,hi in initial])[:,None]
        domain = [endpoints(c,r) for c,r in zip(execution['input_centers'],execution['input_radii'])]
        self._input_lower = np.asarray([[lo for lo,hi in row] for row in domain]).reshape(self._horizon,self._m,1)
        self._input_upper = np.asarray([[hi for lo,hi in row] for row in domain]).reshape(self._horizon,self._m,1)

    @property
    def packet(self):
        return copy.deepcopy(self._packet)

    @property
    def build_info(self):
        return {'source_sha256':self._lib.source_sha256,'compiler_command':self._lib.compiler_command,
                'compiler_version':self._lib.compiler_version,'contract':CONTRACT,
                'platform':sys.platform,'machine':platform.machine(),
                'python':platform.python_version()}

    def simulate(self,initial,inputs):
        import numpy as np
        x = np.ascontiguousarray(initial,dtype=np.float64)
        u = np.ascontiguousarray(inputs,dtype=np.float64)
        if x.ndim == 1:
            x = x[:,None]
        if u.ndim != 3 or u.shape[:2] != (self._horizon,self._m) or u.shape[2] < 1:
            raise ValueError('Wrong input shape or certified horizon')
        count = u.shape[2]
        if x.ndim != 2 or x.shape[0] != self._n or x.shape[1] not in (1,count):
            raise ValueError('Wrong initial shape')
        if not np.isfinite(x).all() or not np.isfinite(u).all():
            raise ValueError('Finite stored inputs required')
        if (x < self._lower).any() or (x > self._upper).any() or (u < self._input_lower).any() or (u > self._input_upper).any():
            raise ValueError('Batch outside certified boxes')
        if x.shape[1] == 1 and count > 1:
            x = np.ascontiguousarray(np.repeat(x,count,axis=1))
        output = np.empty((self._horizon+1,self._p,count),dtype=np.float64)
        pointer = ctypes.POINTER(ctypes.c_double)
        arrays = [v.ctypes.data_as(pointer) for v in (self._e,self._g,self._h,self._d,x,u,output)]
        status = self._lib.run(self._n,self._r,self._m,self._p,self._horizon,count,*arrays)
        if status:
            raise RuntimeError(f'Native execution refused unsupported floating-point environment (code {status})')
        if not np.isfinite(output).all():
            raise OverflowError('Nonfinite native output')
        return output
