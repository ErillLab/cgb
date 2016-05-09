"""
The module for logger initialization.

The logger is configured as specified in configuration file and used across the
package. See Python logging module documentation
(https://docs.python.org/2/library/logging.html) for more information.
"""

import logging
from logging.config import fileConfig

fileConfig('logging.conf')

my_logger = logging.getLogger()
