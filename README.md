# ECG Dilineater and Compressor

This is the course project for the course EE338 - Digital Signal Processing that I did with my friends in a group. In the course mine was Group 1. The other members are [Shlok](https://shlokvaibhav.github.io/) and [Amol](https://www.linkedin.com/in/amol-g-shah/).

---
## Abstract
The analysis of the ECG is a prominent problem in medical science, Most of the useful information
in the ECG is found in the intervals and amplitudes defined by its significant points as shown here:
![ECG](./Dilineation_Schematics.png)

The development of accurate and robust methods
for automatic ECG delineation is a subject of major
importance, especially for the analysis of long
recordings. As a matter of fact, QRS detection
is necessary to determine the heart rate, and as
reference for beat alignment. ECG wave delineation
provides fundamental features (amplitudes
and intervals) to be used in subsequent automatic
analysis.[2]

We use wavelets and filter-banks to perform the
required dilineation as proposed in [2] and [1] to
determine the onset and end of P, T and QRS intervals and determination of R-peak. We also demonstrate
the compression of ECG waves for efficient storage as suggested in [3].
In this assignment we have mainly explored and implemented three research papers -   
[1] “Detection of
ECG characteristic points using wavelet transforms”  
[2] “A Wavelet-Based ECG Delineator: Evaluation on Standard Databases”.  
[3] “ECG Data Compression Using Truncated Singular Value Decomposition”.

More details can be found in the [report](./G1_report.pdf). Also check out the [presentation](./G1_presentation.pptx) if you like.

---
## Code

svd_ECG.m contains the MATLAB code for the compression of ECG signal using SVD. 

The QRS complex detection files are in .m form and .mlx form as well since latter helps in interactive visualization.
The 360 Hz files operate on mitdb and can do only R peak detection and not dilineation since annnotations are only for heart beat
Due to large size of mitdb files, algorithm can take long times to run

The T and P wave detection is available in P_and_T_delineation_250Hz.mlx file. You can run this file independent of others. This file will give the % predictability for P and T wave peak detection.

We have added a tutorial on working with physionet which shows how to import databases into matlab and work with them.
