# galCE
Galaxy chemical evolution code

This script is written to predict the evolution and distribution of elemental abundances in gas, individual stars and synthetic stellar populations in the Milky Way and galaxies based on [Lian et al. 2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474.1143L/abstract), [Lian et al. 202a](https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.2561L/abstract), [Lian et al. 2020b](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.2371L/abstract) and [Lian et al. 2020c](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3557L/abstract). It currently includes sources of metals from AGB stars, Core-collapse supernovae, and type-Ia supernovae. Neutron star mergers and magneto-rotating supernovae are to be included. At the moment, a multi-phase gas accretion and star formation histories are adopted based on Lian et al. 2020b and Lian et al. 2020c. We are working on the code to include more generic (e.g., exponential declining) star formation histories. Please stay tuned for updates! 

To run the GCE model, in the terminal type

$python3 run_mw.py

To check the output of the GCE model, including the model evolution track and simulated mock catalog, run example_multi_burst.ipynb or example_exp.ipynb for examples of different star formation histories. 

If you have any questions, feel free to reach out to me (732592593@qq.com). 

In case this code is used for your study, please cite [Lian et al. 2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474.1143L/abstract) in case of extra-galactic work or [Lian et al. 2020b](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.2371L/abstract) for Galactic work. 
