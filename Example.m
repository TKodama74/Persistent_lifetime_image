%%Example.m
addpath(genpath('programs'))
%%Load samples of integer image and binary mask
load('sample.mat');

qbit = 8; 
sigma = 1;

%%PLT image for connected components (k=0)
k = 0; 
[PLTImg_cc, PLTMap_cc, CentPlot_cc] = Cal_PLTImage(img, roi, qbit, k, sigma);

%%PLT image for hole components (k=1)
qbit = 8; k = 1; sigma = 1;
[PLTImg_hc, PLTMap_hc, CentPlot_hc] = Cal_PLTImage(img, roi, qbit, k, sigma);

%%Display
fig = figure('visible','off');
tiledlayout(2,3);

nexttile;
VOXview(PLTImg_cc,'alpha_def','diff','threshold',1e-5,'colormap',jet);
colorbar;
view([-45 20]);
xlabel('x');
ylabel('y');
zlabel('t');
title("PLT image (CC)");

nexttile;
VOXview(PLTMap_cc,'alpha_def','diff','threshold',1e-5,'colormap',jet);
colorbar;
view([-45 20]);
title("PLT map (CC)");

nexttile;
VOXview(CentPlot_cc,'colormap',jet);
view([-45 20]);
title("Centroid plot (CC)");

nexttile;
VOXview(PLTImg_hc,'alpha_def','diff','threshold',1e-5,'colormap',jet);
colorbar;
view([-45 20]);
title("PLT image (HC)");

nexttile;
VOXview(PLTMap_hc,'alpha_def','diff','threshold',1e-5,'colormap',jet);
colorbar;
view([-45 20]);
title("PLT map (HC)");

nexttile;
VOXview(CentPlot_hc,'colormap',jet);
view([-45 20]);
title("Centroid plot (HC)");

print('-dpng', fullfile('PLTImages.png'));
close(fig);

