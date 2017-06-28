function [EPWu_half, EPDu,EPWi_half, EPDi, hL, PDH, PDW,TP_wi,Tp_w,SCD,SCW,lam_l, IVD] = parameter;
% check the .excel file...

% parameters and correlations from papers
% IVD : intervertebral disc height in order from the first (L1) to the last
% (L5)


% patient height
h = 180 ; %(cm)
age = 36; % age
% According to Taylor evaluate the spine height
CTLh_tot = (h - 47.258)/2.069 ;  % values taken from that paper 


% total IVD height
IVDh = CTLh_tot*0.3 ; 

CTLh_v = CTLh_tot - IVDh ;
% according to Tibbet
hL(1) = 5.55/100 * CTLh_v;
hL(2) = 5.6/100 * CTLh_v;
hL(3) = 5.66/100 * CTLh_v;
hL(4) = 5.63/100 * CTLh_v;
hL(5) = 5.71/100 * CTLh_v;

hL = hL*10; %in mm
%plot(hL, '-.*b'),grid on;
%% EPWu constants
EPWu_a = 1.6836 ;
EPWu_b = -1.598  ;

% you are going to work in axial simmetry, export EPWu/2
EPWu_half = ((EPWu_a.* hL + EPWu_b)./2) ; %it is mm!

%% EPDu constants
EPDu_a = 1.2326;
EPDu_b = 2.838;

EPDu = (EPDu_a .* hL + EPDu_b);

%% EPWi constants
EPWi_a = 1.7618 ;
EPWi_b = -0.7654 ;

% you are going to work in axial simmetry, export EPWu/2
EPWi_half = ((EPWi_a.* hL + EPWi_b)./2) ; %it is mm!

%% EPDi constants
EPDi_a = 1.1348;
EPDi_b = 5.3909;

EPDi = (EPDi_a .* hL + EPDi_b);

%% inferior process width
Ifwi_a = -0.2373;
Ifwi_b = 34.805;

TP_wi = Ifwi_a*hL + Ifwi_b ;

%% transverse process width
Tp_wa = 1.4069;
Tp_wb = 36.851;

%Tp_w = hL*Tp_wa + Tp_wb;
Tp_w = [74.21 88.84 96.89 96.79 97.45];

%% pedicles data
% PDH constants
PDH_a = 0.5528;
PDH_b = 2.0493;

PDH = PDH_a .* hL + PDH_b;
% PDW constants
PDW_a = 0.3683;
PDW_b = 1.2182;

PDW = PDW_a .* hL + PDW_b;
%% spinal canal
SCDa = 0.1216;
SCDb = 17.777;
SCD = SCDa*hL+SCDb;
SCD =SCD*0.9;

SCWa = 0.0905;
SCWb = 16.553;
SCW = SCWa*hL+SCWb;
SCW = SCW*0.9;

lam_l = [36.06 38.58 42.629 36.72 31.766]; %32.382

%% IVD
% IVD according to Shao those are the data for men; for women there are
% other equations
% men
IVD_T12_L1 = 0.519 + 0.004903*age;  %_T12_L1 Actually I won t use it
IVD(1) = 0.680 + 0.006201*age; %_L1_L2
IVD(2) = 0.832 + 0.006687*age; %_L2_L3
IVD(3) = 1.105 + 0.005455*age; %_L3_L4
IVD(4) = 1.076 + 0.006952*age; %_L5_L5
IVD(5) = 0.973 + 0.008630*age; %_L5_S1

IVD = IVD *10;%in mm
%{
IDV2(1)= 12.35/100*IVDh  ;%L1-L2
IDV2(2)= 15.07/100*IVDh   ;%L2-L3
IDV2(3)= 18.41/100*IVDh   ;%L3-L4
IDV2(4)= 18.83/100*IVDh   ;%L4-L5
%IDV2(5)=
%}
% women
%{
IVD_T12_L1 = 0.433 + 0.004903*age;  %_T12_L1 Actually I won t use it
IVD(1) = 0.627 + 0.004840*age; %_L1_L2
IVD(2) = 0.817 + 0.004771*age; %_L2_L3
IVD(3) = 0.985 + 0.004982*age; %_L3_L4
IVD(4) = 1.051 + 0.005052*age; %_L5_L5
IVD(5) = 0.926 + 0.008170*age; %_L5_S1
%}







