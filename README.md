# ColleciveDynamicKModule
On the scale entwined dynaimcs emerged in the potassium channel module.
======================================================================
This is the codes and data fro the paper titled: "On the scale-entwined kinetics emerged in the voltage-gated potassium channel module".


Directory  "HybridGelspStepbyStep"
Here put the program "PowerSpect_par.m", to calculate the power spectra of the hybrid process. It is the microscopic model. Channel state transition is discrete, jump.
Voltege V is continous.

The advanced Giellespie algrithm is used.
"FindLoca2.m" is the function that will be called in the main program.

"Samplepath.m"  , to calculate the sample path of the hybrid process.

In the directory "Data", put  the calculated data to plot the figures in the paper.
"plot(Freq,cumsum(AV_df_V))"  to plot the Cumulative power spectrum;  " plot(Freq,AV_df_V)", to plot the power spectral density.

All the parameer values are saved in these ".mat"  data files.

+==============================================

Directory: "DiffusThrePart"
Here put the program "PowerSpect_par", to calcuLate the power spectra of the diffusion appximation model. The algrithm is Eular and Mote Calor.

"Samplepath.m": calculate the sample paths of the diffusion system.
In the directory "Data", put  the calculated data to plot the figures in the paper.
"plot(Freq,cumsum(AV_df_V))"  to plot the Cumulative power spectrum;  " plot(Freq,AV_df_V)", to plot the power spectral density.

All the parameer values are saved in these ".mat"  data files.
