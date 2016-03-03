import logging
from logging.config import fileConfig

fileConfig('logging.conf')

my_logger = logging.getLogger()
