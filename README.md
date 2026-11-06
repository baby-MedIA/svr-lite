SVR-lite - Slice-to-volume reconstruction software for fetal MRI
====================
SVR-lite is a simplified, refactored and optimised software package for 3D and gated 3D+T/4D rigid slice-to-volume reconstruction.


Installation
------------

Installation instructions are available in InstallationInstructions.txt

SVR-lite algorithms
--------------------

SVR-lite implements rigid slice-to-volume reconstruction, based on the algorithm described here:
> Kuklisova-Murgasova, M., Quaghebeur, G., Rutherford, M. A., Hajnal, J. V., & Schnabel, J. A. (2012). Reconstruction of fetal brain MRI with intensity matching and complete outlier removal. Medical Image Analysis, 16(8), 1550–1564.: https://doi.org/10.1016/j.media.2012.07.004

SVR-lite also implements an extension of the original 3D algorithm to 3D+time/4D data described here:
> Boutillon, A., Clarke, N., Woodgate, T., West, D., Schneider, A., Franklin, R., Price, A., Hajnal, J., Pushparajah, K., Lloyd, D., & Deprez, M. (2026). SPARC: Slice-to-volume Pipeline for Automated Reconstruction of gated 3D+time fetal Cardiac MRI. arXiv preprint arXiv:2608.18616. https://doi.org/10.48550/arXiv.2608.18616

The 3D+time/4D reconstruction algorithm was inspired by 4D cardiac magnitude reconstruction described here:
> van Amerom, J. F. P., Lloyd, D. F. A., Deprez, M., Price, A. N., Malik, S. J., Pushparajah, K., van Poppel, M. P. M, Rutherford, M. A., Razavi, R., Hajnal, J. V. (2019). Fetal whole-heart 4D imaging using motion-corrected multi-planar real-time MRI. Magnetic Resonance in Medicine, 82(3): 1055-1072. : https://doi.org/10.1002/mrm.27858

SVR-lite implementation
------------------------
The reconstruction functions in SVR-lite are derived from the original 3D rigid implementation in https://gitlab.com/mariadeprez/irtk-simple, where image handling and registration functions were derived from IRTK (https://github.com/BioMedIA/IRTK).

3D+T/4D rigid implementation was inspired by https://github.com/jfpva/fetal_cmr_4d.

SVR-lite is derived from a subset of SVRTK code at https://github.com/SVRTK/SVRTK/tree/fbd77ab0a1d4febc47dda349a06283b3c3589754. 

SVR-lite uses image handling and registration functions from MIRTK ( https://github.com/BioMedIA/MIRTK), a successor of the original IRTK (https://github.com/BioMedIA/IRTK). MIRTK was originally adapted for SVRTK in the SVRTK repository https://github.com/SVRTK/SVRTK.
 
SVR-lite includes a copy of MIRTK based on the version adapted for SVRTK (https://github.com/SVRTK/MIRTK).
