# CARMER
Carbon-Mercury ocean-atmosphere model

As used in Dal Corso et al. 2020 Nat. Comm.

Contents:
CARMER_frontend.m // Model frontend file // run this code to solve model and plot results // sets parameters and starting values, runs solver
CARMER_equations.m // Model equations file // do not run this code directly // contains flux and reservoir calculations
d13c_hg_data.mat // Datafile containing C and Hg data // plots automatically against model

Reproduction of paper plots:
Lines 40-140 in the equations file cover the three scenarios tested in the paper. Comment each code block in our out to run.

Notes:
Requires MATLAB. Tested in R2018a on win10 x64. Run time << 1s. No installation required.

Fixes:
This version includes a bug fix  - the isotopic composition of carbon removed from the atmosphere via carbonate weathering was taken to be an average of the carbonate rock and atmospheric carbon, rather than being atmospheric carbon alone. This acts to change the steady state C isotope value slightly, requiring a slightly different choice of the assumed terrestrial productivity to reproduce the paper plots exactly. The paper conclusions are unaffected.
