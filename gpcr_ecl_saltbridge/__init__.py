"""
gpcr_ecl_saltbridge
Identify GPCR structures with saltbridges between extracellular loops 2 and 3.
"""

# Add imports here
from .gpcr_ecl_saltbridge import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
