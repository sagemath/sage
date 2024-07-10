from sage.misc.lazy_import import lazy_import

lazy_import("sage.calculus.transforms.fft", ["FastFourierTransform", "FFT"])
lazy_import("sage.calculus.transforms.dwt", ["WaveletTransform", "DWT"])
from sage.calculus.transforms.dft import IndexedSequence
del lazy_import
