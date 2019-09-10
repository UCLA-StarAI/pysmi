
import logging
from sys import stdout

logger = logging.getLogger("pysmi")
#logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler(stdout)
logger.addHandler(stream_handler)
