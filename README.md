# ABSTRACT #
ECG signal conditioning by morphological Filtering

## MATLAB (v.R2022a) Requirements ## 
- **Communications Toolbox** 
- **Image Processing Toolbox** (*if you want to use the build-in functions like: imopen, imclose, imerode and imdilate to process the signal*)

## STRUCTURE ##
- **ECG_MMF**: main class of the project.
- **Gen_Drift**, **Gen_Noise** and **Gen_Strel**: generation of structures and artifacts useful for the froject.
- **erosion.m**, **dilatation.m**, **closing.m** and **opening.m**: definition of the morphological operators.

## INFO ##
The directory '*dataset*' contains a series of ECG signals recovered from the physionet database (LINK: https://physionet.org/content/mitdb/1.0.0/)

The functions '*erosion.m*', '*dilatation.m*', '*closing.m*' and '*opening.m*' are used to process the signals. The functions '*GenDrift.m*', '*GenStrel.m*' are used to generate respectively a baseline drift in the signal and to generate a bunch of Structuring Element to operate on signals. 

If you do not want to use the **Communication Toolbox** with the '*awgn*' (Add White Gaussian Noise) function, instead you can adopt the '*GenNoise*' function to generate noise.

## DOCUMENTATION ##
The two pdf file, *'paper4.pdf'* and *'extra_for_paper4.pdf'*, are the papers used to elaborate the project. The last pdf file named as Progect_Review consist in a simple presentation of the project.
