# SIDM_dark_star_code
Code for dark stars powered by Self-Interacting Dark Matter (SIDM). 
In case nucrat.so does not work, run the following command in the shell to rebuild the nucrat package from the fortran script with f2py:
f2py -c --fcompiler='gfortran' -m nucrat nucrat.f
A module with a name like 'nucrat.cpython-38-darwin.so' will be created. Rename this module into 'nucrat.so' and you are all set!
