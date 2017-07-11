function [EPWu_half, EPDu, EPWi, EPDi, hL, PDH, PDW, IVD, lam_l, sc_d, sc_w, TP_wu, TP_wi, F_h] = parameter_visual_human;
% I measured the dimensions from the model reconstracted at the University
% of Melbourne
%check the acronyms in the file
%C:\Users\lavecchc\Dropbox\Internship_Melbourne\dimensions_vh.xls

EPWu = [46.598 49.738 48.967 52.687 54.498]; % there is also the T12= 44.898 
EPWu_half = EPWu/2;

EPDu = [36.159 43.51 42.527 40.23 37.941]; %T12= 37.4686 

PDH = [ 15.719 16.278 17.416 18.045 17.256 ] ;% 19.2016
PDW = [ 8.189 9.288 13.597 14.963 15 ]; % 8.821

hL = [ 29.404 31.625 30.719 30.129 28.976]; %28.613 

IVD = [9.065 11.112 12.911 11.072 9.279 ]; %

lam_l = [36.06 38.58 42.629 36.72 31.766]; %32.382

sc_d = [18 16.903 16.035 17 17];%19.256 %or 19.644 20.436 

sc_w = [22.265 22.244 22.564 24 25]; %20.398 %orig 26.043 29.743

TP_wu = [74.214 88.84 96.893 96.785 97.451]; %51.568

TP_wi = [24.97 22.429 26.8 31.254 42.135];%23.039

EPWi = [ 48.84 50.47 52.186 53.437 54.1929]; %46.959

EPDi = [ 37.72 44.754 39.037 38.729 36.127]; %35.681

F_h = [47.954 53.887 53.793 56.049 48.439]; %45.441 %facet length from the upper surface of the vertebral body

end