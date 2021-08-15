# Extract power line from mobile LiDAR point cloud
The mobile laser scanning (MLS) system can quickly collect the point cloud data of power lines and power line corridors that lay along urban road. Accurate and efficient extraction of power lines from the point cloud is the basis of power line inspections and maintenance. This project presents a complete algorithm for power line extraction and modelling from MLS point clouds. The algorithm first extracted candidate power line points from non-ground points based on the analysis of linear feature. The catenary model is used to model and densify power lines. The sample data provided by the International Society for Photogrammetry and Remote Sensing Working Group (ISPRS WG) III/5 was used to test the performance of the method.

Contents
****

-   [Compile](#Compile)
-   [Examples](#Examples)
-   [Our publications](#our-publications)
-   [Results](#Results)

Compile
--------------
Compile in the matlab with:
```cpp
%build from source

mex -I'.\Eigen\eigen3' extractPLs.cpp detectpowerline.cpp kdtree.cpp 3dKDtree.cpp eigenmatrix.cpp
```
Example
--------------
Run in Matlab: demo_extract_powerline.m
key codes:
```cpp
radius = 0.5;
angleThr = 10;
LThr = 0.98;
tic
[isPLIndex] = extractPLs(nonGroundPoints,radius,angleThr,LThr);
toc
isPLIndex = logical(isPLIndex);

figure
pcshow(nonGroundPoints(isPLIndex,:))
title('candidate powerline points')
view(60,25)
```

Results
--------------
<div align=center>
<img src="https://github.com/zwshi-pku/3DLiDAR/blob/main/ExtractPowerLine/re_f1_nonGroundPoints.png"  height="70%" width="70%">
<img src="https://github.com/zwshi-pku/3DLiDAR/blob/main/ExtractPowerLine/re_f2_candidate%20powerline%20points.png"  height="70%" width="70%">
<img src="https://github.com/zwshi-pku/3DLiDAR/blob/main/ExtractPowerLine/re_f3_candidate%20powerline%20points%20clusters.png"  height="70%" width="70%">
<img src="https://github.com/zwshi-pku/3DLiDAR/blob/main/ExtractPowerLine/re_f4_powerline%20points%20clusters.png"  height="70%" width="70%">
<img src="https://github.com/zwshi-pku/3DLiDAR/blob/main/ExtractPowerLine/re_f5_colorization%20clusters.png"  height="70%" width="70%">
<img src="https://github.com/zwshi-pku/3DLiDAR/blob/main/ExtractPowerLine/re_f6_powerLines%20clusters.png"  height="70%" width="70%">
<img src="https://github.com/zwshi-pku/3DLiDAR/blob/main/ExtractPowerLine/re_f7_Power%20line%20model.png"  height="70%" width="70%">
</div>


Our Publications
--------------

Please consider citing our papers if you find these codes help your research.

**Zhenwei Shi**, Yi Lin, and Hui Li. **Extraction of urban power lines and potential hazard analysis from mobile laser scanning point clouds.** International Journal of Remote Sensing 41, no. 9 (2020): 3411-3428.

**Zhenwei Shi**, Zhizhong Kang, Yi Lin, Yu Liu, and Wei Chen. **Automatic recognition of pole-like objects from mobile laser scanning point clouds.** Remote Sensing 10, no. 12 (2018): 1891.



