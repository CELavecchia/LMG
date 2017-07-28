function [dimensions] = averageModelInput(varargin)


% check the .excel file...

% parameters and correlations from papers
% IVD : intervertebral disc height in order from the first (L1) to the last
% (L5)

%% parse input
h = varargin{1}.height; % patient height (cm)
age = varargin{1}.age; % Patient age (years)
sex = varargin{1}.sex;

% this one can be used in case of hybrid models. Presently, all the
%  dimensions are required from the scan.

% switch varargin{1}.type
%     case NaN        
%         h = 180; % patient height (cm)
%         age = 36; % Patient age (years)
%         sex = 'm';
%     case 0
%         h = varargin{1}.height; % patient height (cm)
%         age = varargin{1}.age; % Patient age (years)
%         sex = varargin{1}.sex;
%  
%     case 1 %subject specific dimensions
%         h = varargin{1}.height;% patient height (cm)
%         age = varargin{1}.age; % Patient age (years)
%         sex = varargin{1}.sex;
%        
% end

%%

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
    
if strcmp(sex,'m')==1
    % men
    IVD_T12_L1 = 0.519 + 0.004903*age;  %_T12_L1 Actually I won t use it
    IVD(1) = 0.680 + 0.006201*age; %_L1_L2
    IVD(2) = 0.832 + 0.006687*age; %_L2_L3
    IVD(3) = 1.105 + 0.005455*age; %_L3_L4
    IVD(4) = 1.076 + 0.006952*age; %_L5_L5
    IVD(5) = 0.973 + 0.008630*age; %_L5_S1
    
    IVD = IVD *10;%in mm
else
    % women

    IVD_T12_L1 = 0.433 + 0.004903*age;  %_T12_L1 Actually I won t use it
    IVD(1) = 0.627 + 0.004840*age; %_L1_L2
    IVD(2) = 0.817 + 0.004771*age; %_L2_L3
    IVD(3) = 0.985 + 0.004982*age; %_L3_L4
    IVD(4) = 1.051 + 0.005052*age; %_L5_L5
    IVD(5) = 0.926 + 0.008170*age; %_L5_S1
    IVD = IVD *10;%in mm
end


%% pedicle angles
%trasverse
PDt_a = 0.042;
PDt_b = 4.683;
PDt = PDt_a .* hL + PDt_b;

%sagittal
PDs_a = -1.246;
PDs_b = 46.075;
PDs = PDs_a .* hL + PDs_b;

%% Create output structure

dimensions.EPWu_half =EPWu_half ;
dimensions.EPDu= EPDu; 
dimensions.EPWi=EPWi_half;
dimensions.EPDi=EPDi;
dimensions.hL=hL; 
dimensions.PDH=PDH; 
dimensions.PDW=PDW;
dimensions.IVD=IVD;
dimensions.lam_l=lam_l; 
dimensions.sc_d=SCD;
dimensions.sc_w=SCW;
dimensions.TP_wu=Tp_w;
dimensions.TP_wi=TP_wi;
dimensions.PDt=PDt;
dimensions.PDs=PDs;

