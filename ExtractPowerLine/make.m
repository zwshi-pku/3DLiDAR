% /****************************************************************************
%  * Please cite the following paper, If you use this code in your work.
%  *
%  * Zhenwei Shi, Yi Lin, and Hui Li. "Extraction of urban power lines and 
%  * potential hazard analysis from mobile laser scanning point clouds." 
%  * International Journal of Remote Sensing 41, no. 9 (2020): 3411-3428.
%  *
%  * The paper can be downloaded from
%  * https://www.tandfonline.com/doi/full/10.1080/01431161.2019.1701726
%  * Copyright
%  * Institute of Remote Sensing & GIS, 
%  * School of Earth and Space Sciences, Peking University (www.pku.edu.cn)
%  * Zhenwei Shi; Yi Lin; Hui Li
%  * contact us: zwshi@pku.edu.cn; lihui@pku.edu.cn
% *****************************************************************************/

%Compile in the matlab with:
%build from source
mex -I'.\Eigen\eigen3' extractPLs.cpp detectpowerline.cpp kdtree.cpp 3dKDtree.cpp eigenmatrix.cpp