"""
The module for logger initialization.

The logger is configured as specified in configuration file and used across the
package. See Python logging module documentation
(https://docs.python.org/2/library/logging.html) for more information.
"""
import os
import logging
from logging.config import fileConfig

fileConfig(os.path.join(os.path.dirname(__file__), 'logging.conf'))

my_logger = logging.getLogger()
