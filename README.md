#ImperialSRUS

This code provides a simplified implementation of the framework in below paper. Deconvolution in the paper is replaced by normalised cross-correlation in this code for localisation.

J. Yan, T. Zhang, J. Broughton-Venner, P. Huang and M. -X. Tang, "Super-Resolution Ultrasound Through Sparsity-Based Deconvolution and Multi-Feature Tracking," in IEEE Transactions on Medical Imaging, vol. 41, no. 8, pp. 1938-1947, Aug. 2022, doi: 10.1109/TMI.2022.3152396.  


Code is published under a Creative Common license for non-commercial use (CC-BY-NC), and therefore can be used for non-commercial, personal or academic use as long as the paper is correctly cited.

Code is developed and tested in Matlab2021a on a desktop (CPU: AMD Ryzen 9 5900 Processor, GPU: Nvidia Geforce RTX3080, RAM: 128 Gb). 
Pixel map of input image should be dense enough to make it able to estimate point spread function (PSF) from the CEUS sequence. According to our experience, the PSF estimation works well when each dimension of PSF is at least 5 pixels. 
