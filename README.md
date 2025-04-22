# GRL-2025
Post-processing scripts for "Flow-induced Compaction Localization in a Rock-like Medium" by Arnold Bachrach and Yaniv Edery - Technion - Israel Institute of Technology

The data files required to reproduce figures 2 and 3 of the main text and figures 1 and 2 of the supplementary material is available here, along with the MATLAB scripts used to post-process it. The MATLAB scripts used for post-processing are available also at: https://github.com/Yaniv-Edery/GRL-2025. The initial PIV analysis was performed using PIVlab 3.09, an open-source MATLAB-based tool accessible at https://github.com/Shrediquette/PIVlab. MATLAB R2022b was used for both the PIVlab and the post-processing scripts. For each experiment ('S1','S2','S3') we provide three datasets which are used for three different post-processing scripts:


* 'Workspace strain stress updated' includes the measured injection pressure ('P new' $Pa$), the measured discharge ('Q new' $cm^3/s$) and the displacements in $x_1$ and $x_2$ ('u original' and 'v original' $m$, respectively ), which where calculated by the PIVlab with a reference frame at $P=0.15$ $MPa$. The displacement fields are provided in a temporal resolution of 0.1 s. This workspace is used with the script named 'Code strain stress' for calculating the averaged strain-stress curve and the strains along $x_1$ (figure 2 'a' and 'b', respectively). Within this script you can also find a section that plots the fluid flux versus $P$ and calculates the maximum Reynolds number throughout the experiment.  

* 'Workspace 5 sec window updated' includes the displacements in $x_1$ and $x_2$ ('u original' and 'v original' $m$, respectively), which where calculated by the PIVlab with a reference frame at $P=0.15 MPa$. In this dataset we provide the displacements in a temporal resolution of 0.01 s. Yet, this dataset include only 5 s of the experiment, around the localization. This workspace is used with the script named 'Code 5 sec window'  for calculating the 5 s window curve in figure 3 'a'.

* 'Workspace pre localization ref updated' includes the displacements in $x_1$ and $x_2$ ('u original' and 'v original' $m$, respectively), which where calculated by the PIVlab with the reference frame being the last frame before the onset of the localization, at a temporal resolution of 0.01 s, and along 3 s of the localization process. This workspace is used with the script named 'Code pre localization ref', for calculating 2D displacement gradient fields (figure 3 'b' and figure 2 in the supplementary material) as well as 2D displacement fields at the onset of localization (figure 1 in the supplementary material).    


Cite as: Bachrach, A., & Edery, Y. (2025a). Flow-induced compaction localization in a rock-like medium.
