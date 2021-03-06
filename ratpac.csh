#!/bin/csh

setenv RATROOT /home/enrg/Software/WMUtils/ratpac/build
setenv RATSHARE /home/enrg/Software/WMUtils/ratpac
setenv PATH "$RATROOT/bin:$PATH"
if ({$?LD_LIBRARY_PATH}) then
  setenv LD_LIBRARY_PATH "${RATROOT}/lib:$LD_LIBRARY_PATH"
else
  setenv LD_LIBRARY_PATH "${RATROOT}/lib"
endif
setenv ROOT_INCLUDE_PATH "$RATROOT/include"

# For Mac OS X
if ({$?DYLD_LIBRARY_PATH}) then
  setenv DYLD_LIBRARY_PATH "${RATROOT}/lib:\$DYLD_LIBRARY_PATH"
else
  setenv DYLD_LIBRARY_PATH "${RATROOT}/lib"
endif

if ({$?PYTHONPATH}) then
  setenv PYTHONPATH "$RATSHARE/python:$PYTHONPATH"
else
  setenv PYTHONPATH "$RATSHARE/python"
endif
setenv GLG4DATA "$RATSHARE/data"
