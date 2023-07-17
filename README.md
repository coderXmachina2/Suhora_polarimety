# Suhora_polarimety data processing pipeline

![XL1n2](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/MISSION_LOGO.png)

R band optical polarimetry of EE Cephei with birefringent Wollaston prism/ birefringent Savart plate polarimetric filter. This code encompasses the tools and processing pipe to compute polarization degree and position angle as a function of time with addititonal light curve visualization capable.

The objective of the experiment is to investigate the variation of Polarization Degree (PD) and Polarization Angle (PA) as a function of time in conjunction with the 2020 Eclipse of the debris disk of the variable star EE Cep (d `≈` 2.75 kpc). The eclipse standardized magnitude is visualized in the light_curve below where the ingress occuring in February and the eggress in May-June.

![light_curve](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/EE_Cep_Stand_Mag_flux_plot.PNG)

Published as a European Astronomical Society (EAS) conference poster in 2022. Link found here:

https://k-poster.kuoni-congress.info/eas2022/poster/1694cb25-448d-4d5e-a953-80b932538b2a

Custom code was written to reduce, process, and analyze data from the 60 cm telescope on Mt. Suhora Poland. The raw data (FITS files) are not included in this repo, however, the stats containing Stokes polarimetry values are included. Derivation of the statistics are as follows with in depth documentation coming soon:

Steps:
1. Reduce CCD image. Subtract the image with bias, subtract the image with dark, divide image by bias subtracted master flat. Corrects the flux distribution of the image. Flat fielding reduces the interpixel variations in the image (the bright spot in the image center).

![CCD_Reduction](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/ccd_reduction.png)

2. Perform photometric counts extraction (using circular aperture photometry) on peak detected sources within an input defined region of interest (target and polarization standards). Aperture size must be adjusted manually.

![Thresholding](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/thresholding.png)

![Region_Of_Interest](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/region_of_interest.png)

![RTP_1](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/Run_the_pipe_1.png)

![RTP_2](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/Run_the_pipe_2.png)

4. Perform computation of normalized Stokes parameters q and u from source counts based on equations from Ramanprakash et al. 2019 where q = (I<sub>0</sub>-I<sub>90</sub>)/(I<sub>0</sub>+I<sub>90</sub>) and u = (I<sub>45</sub>-I<sub>135</sub>)/(I<sub>45</sub>+I<sub>135</sub>). As beautiful as the images are they all end up as a table full of numbers. And that is where the true beauty lies.

![XL1n2](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/Excels_1_and_2.PNG)

![XL1n2_combined](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/Excels_1_and_2_combined.PNG)

6. Visualize Stokes parameters. Full plot and single plot options are available. Full scatter and mean scatter options are available.

![master_pol](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/master_plot_polarimetry.png)

![full_q_u](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/target_full_q_u.png)

![mean_q_u](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/target_mean_q_u.png)

8. Check for polarimetric stability/ time variability of Stokes parameters and calculate instrumental polarization

![q_stab](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/q_stability.png)

![u_stab](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/u_stability.png)

10. Calculate Polarization Degree (PD) and Polarization Angle (PA)

![PD_stab](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/PD_stab.png)

![PA_stab](https://github.com/coderXmachina2/Suhora_polarimety/blob/main/github_imgs/PA_stab.png)


