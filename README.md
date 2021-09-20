# Suhora_polarimety
Suhora Polarimetry

R band optical polarimetry of EE Cephei with birefringent Wollaston prism/ birefringent Savart plate polarimetric filter. Basically the tools and processing pipe to compute polarization degree and position angle

Steps:
1. Reduce CCD image. Subtract bias, subtract dark, divide by bias subtracted flat.
2. Perform photometry counts extraction on peak detected source
3. Perform Polarimetric computation on 4 source counts using equations from Ramanprakash et al. 2019.

Data not included in Repo
