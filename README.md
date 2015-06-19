# my5spec
Utility to generate dynamic spectra

Usage: ``my5spec [-a aver_time] [-n nchan] [-l time_limit] [-o offset] INFILE FORMAT OUTFILE``

``INFILE``    - the name of the input file

``FORMAT``    - mark5access data format in form ``<FORMAT>-<Mbps>-<nchan>-<nbit>``

``OUTFILE``   - basename for the output files. Output files will be called `OUTFILE_n`, where n is channel number

Parameters:

  ``-a aver_time``  - approximate integration time per spectrum in milliseconds (1 ms)
  
  ``-n nchan``      - number of spectral channels (128)
  
  ``-l time_limit`` - total time in seconds
  
  ``-o offset``     - offset in seconds from the file beginnig (0)
