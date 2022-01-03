# Suhora_polarimety data processing pipeline

R band optical polarimetry of EE Cephei with birefringent Wollaston prism/ birefringent Savart plate polarimetric filter. Basically the tools and processing pipe to compute polarization degree and position angle

Steps:
1. Reduce CCD image. Subtract bias, subtract dark, divide by bias subtracted flat.
2. Perform photometric counts extraction on peak detected sources (target and polarization standards)
3. Perform computation normalized Stokes parameters from 4 source counts based on equations from Ramanprakash et al. 2019.
4. Visualize Stoles parameters
5. Calibrate for instrumental polarization
6. Calculate Polarization Degree and Polarization Angle (PA)

Excel files included in repo
