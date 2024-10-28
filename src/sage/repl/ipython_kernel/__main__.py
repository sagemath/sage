# sage_setup: distribution = sagemath-repl
from ipykernel.kernelapp import IPKernelApp
from sage.repl.ipython_kernel.kernel import SageKernel
IPKernelApp.launch_instance(kernel_class=SageKernel)
