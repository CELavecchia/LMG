clc; clear all; close all;
%function to call the gui and start the model
%w = warning ('off','all');
%rmpath('./output')
%warning(w)

%[model, var] =GUI_test1
%{
get:
- type of model (a.average, b.subject-specific)
- num of bodies (you can use the whole spine or only a functional unit.. it
has to be improved)
- output (.stl,.feb,.inp)
- alpha (lumbar curvature angle, if zero then the bodies are arranged in vertical)

%}

nbodies =2;

%  [dimensions] = parameter_average(model.average.age,model.average.height
%  ,model.average.sex); % from output gui

[EPWu_half, EPDu,EPWi_half, EPDi, hL, PDH, PDW,TP_wi,Tp_w,SCD,SCW,lam_l, IVD] = parameter_no_input;

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

% fitting and parameterize the model
mesh_struct_IVD2 = model_elaboration_gui(dimensions);

% Mesh vertebrae
%nbodies,stl
stl = 0;
Lmes = FU_whole_model_gui(2,stl);

%orientation post mesh
orientation_gui(Lmes,nbodies)%, alpha); % check it!

plot3(Lmes(1).VT(:,1),Lmes(1).VT(:,2),Lmes(1).VT(:,3),'.r'),hold on;

plot3(Lmes(2).VT(:,1),Lmes(2).VT(:,2),Lmes(2).VT(:,3),'.b'),hold on;

% write FEBio
%name = ('./output/feb/model_febio');
%FEB_struct = febio_preprocessing_gui(name,nbodies, Lmes,mesh_struct_IVD2)
