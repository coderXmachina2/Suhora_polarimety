# Suhora_polarimety data processing pipeline

R band optical polarimetry of EE Cephei with birefringent Wollaston prism/ birefringent Savart plate polarimetric filter. Basically the tools and processing pipe to compute polarization degree and position angle

Steps:
1. Reduce CCD image. Subtract bias, subtract dark, divide by bias subtracted flat. Corrects the flux distribution of the image. Reduces interpixel variations in the image (the bright spot in the middle).

![CCD_Reduction](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/ccd_reduction.png)

2. Perform photometric counts extraction on peak detected sources within an input defined region of interest (target and polarization standards).

![Thresholding](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/thresholding.png)

![Region_Of_Interest](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/region_of_interest.png)

4. Perform computation of normalized Stokes parameters q and u from source counts based on equations from Ramanprakash et al. 2019 where q = ()/() and u = ()/()

![XL1n2](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/Excels_1_and_2.PNG)

![XL1n2_combined](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/Excels_1_and_2_combined.PNG)

6. Visualize Stokes parameters
7. Calibrate for instrumental polarization
8. Calculate Polarization Degree and Polarization Angle (PA)

Excel files included in repo
