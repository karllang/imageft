# imageft

C. Glotzbach (07.03.2017)

How to get the FT-correction from microscopic images taken parallel and perpendicular to the crystallographic c-axis of apatite or zircon grains:

1. Open the script batch_Ft_grain.m in the Matlab Editor
2. Change variables (line 3-13), e.g. adjust the resolution of images (microns per pixel) and if available change isotope ratios (r232_238 and r147_238) to those measured. If not available leave zero (average values will be used to estimate the Ft). Note that the final Ft-value may differ from this estimate, if isotope ratios differ from the average values.
2. Run the script 
3. It will automatically trace the grain outline, you may have to adjust it
4. The orientation of the crystallographic c-axis will be shown, you may have to adjust it by moving the end-points. Double click the one of the end-points if you are happy.
5. The outline of the grain perpendicular to the crystallographic c-axis will be traced. Adjust it if necessary.
6. The Ft-values of broken grains can be estimated. If this is applicable to your grain, you can use the mouse to add some missing pieces of the grain (draw a polygon overlapping with the grain outline). Looks like this:

Note that the FT is calculated assuming some mean U/Th/Sm ratios. This is good to get an idea of the FT. To get the correct FT required for the final age calculation you have to change ‚r232_238‘ and ‚r147_238‘ to measured values derived from the ICP-MS.

Led me know if you have any idea to make this code more user-friendly. 
