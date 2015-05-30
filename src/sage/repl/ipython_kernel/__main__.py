from IPython.kernel.zmq.kernelapp import IPKernelApp
from sage.repl.ipython_kernel.kernel import SageKernel
IPKernelApp.launch_instance(kernel_class=SageKernel)
