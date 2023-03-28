# Wavenumber-calibration
The code and data relate to the following paper:
Dongyue Liu, Bryan M. Hennelly, "Wavenumber Calibration by Modelling the Raman Spectrometer" (under review JRS 2023)
See paper for further details.

 

The code provided is just one example for 4-acetamidophenol material at 300 lines/mm.

 

The user should first download the code and extract the data from the folder "reference wavenumber spectra" folder such that it can be accessed by the matlab files. The various m files are explained as follows:

 

Run ace_g300_leave_all_out.m to see the results of Leave all out method. This code inlcudes the code for the algorithm described in the paper, with the various steps being clearly defined.

 

Run ace_g300_leave_one_out.m to see the results of Leave one out method. This code inlcudes the code for the algorithm described in the paper, with the various steps being clearly defined.

 

Run ace_g300_leave_left_out.m and ace_g300_leave_right_out.m can be called to see the results for Leave half out method. This code inlcudes the code for the algorithm described in the paper, with the various steps being clearly defined.

 

fpeak.m is the code used to find the peak positions in the spectrum
 

lorentzfit.m is the code used to fit a lorenzian to a found peak so that the position can be determined with hsub-pixel accuracy


rawdata.mat contains the data (100 spectra)
