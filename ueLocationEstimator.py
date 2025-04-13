####################################################################################################################################################################
# Importing Libraries
import pip

try:
    import os
except:
    pip.main(['install', '--user', os])
    import os

def install_or_import(packageList):
    for package in packageList:
        try:
            _import_(package)
        except:
            os.system("pip3 install " + package)

    return

packageList = ["numpy"]
install_or_import(packageList)

import numpy as np
####################################################################################################################################################################

# Importing System Config
sysConfig = getSystemConfig()

# Fetching 