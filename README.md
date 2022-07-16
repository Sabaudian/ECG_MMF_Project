### ABSTRACT ###
ECG signal conditioning by morphological Filtering

### MATLAB_R2022a Requirements ### 
Communications Toolbox (Image Processing Toolbox if you want to use the buid-in functions like: imopen, imclose, imerode and imdilate to process the signal) 

The directory 'dataset' containing a series of ECG signals recovered from the physionet database (LINK: https://physionet.org/content/mitdb/1.0.0/)

The functions 'erosion.m', 'dilatation.m', 'closing.m' and 'opening.m' are used to process the signals. The functions 'GenDrift', 'GenStrel' are used to generate respectively a baseline drift in the signal and to generate a bunch of Structuring Element to operate on signals. 

The two pdf file are the documentations used to elaborate the project.
