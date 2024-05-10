from sage.misc.lazy_import import lazy_import

lazy_import("sage.calculus.transforms.fft", ["FastFourierTransform", "FFT"])
lazy_import("sage.calculus.transforms.dwt", ["WaveletTransform", "DWT"])
from .dft import IndexedSequence
