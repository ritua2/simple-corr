|Travis Build|

simple-corr
##############

A repo containing a simple correlation function code (and tests)

To compile the CUDA version of the code:

nvcc -arch=sm_35 -o simple-corr-version3 simple-corr-version3.cu  

.. |Travis Build| image:: https://travis-ci.com/manodeep/simple-corr.svg?branch=master
   :target: https://travis-ci.com/manodeep/simple-corr
   :alt: Build Status
