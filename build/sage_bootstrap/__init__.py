
from sage_bootstrap.config import Configuration
config = Configuration()

from sage_bootstrap.stdio import init_streams
init_streams(config)

from sage_bootstrap.logger import init_logger
init_logger(config)
