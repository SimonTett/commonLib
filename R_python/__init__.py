import rpy2.situation
import os
import logging
os.environ['R_HOME'] = os.path.expandvars(rpy2.situation.get_r_home())  # from https://github.com/rpy2/rpy2/issues/796
logging.info(f"Set R_HOME to {os.environ['R_HOME']}")