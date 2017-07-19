clc; clear all; close all;
fontSize=15;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
%function to call the gui and start the model
fprintf('--------------------LMG: LUMBAR MODEL GENERATOR------------------\n\n\n')

% file name
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
outputpathName=fullfile(defaultFolder,'data','output','feb');
febName=fullfile(outputpathName,'model_febio');


%[model, var] =GUI_test1
%{
get:
- type of model (a.average, b.subject-specific)
- num of bodies (you can use the whole spine or only a functional unit.. it
has to be improved)
- output (.stl,.feb,.inp)
- alpha (lumbar curvature angle, if zero then the bodies are arranged in vertical)

%}

nbodies =5;

fprintf('-----------------Build the geometrical model--------------------\n');
%  [dimensions] = parameter_average(model.average.age,model.average.height
%  ,model.average.sex); % from output gui

[EPWu_half, EPDu,EPWi_half, EPDi, hL, PDH, PDW,TP_wi,Tp_w,SCD,SCW,lam_l, IVD, PDt, PDs] = parameter_no_input; %add output

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

alpha = 43.49;

% fitting and parameterize the model
mesh_struct_IVD2 = model_elaboration_gui(dimensions);

cFigure;
%subplot(1,2,1);

 title('Solid tetrahedral meshing IVD','FontSize',fontSize);
 xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
 patch('Faces',mesh_struct_IVD2(1).Fb,'Vertices',mesh_struct_IVD2(1).V,'FaceColor','flat','CData',mesh_struct_IVD2(1).faceBoundaryMarker,'lineWidth',0.2,'edgeColor',edgeColor);
axis tight;  axis equal;  %grid on; 
camlight headlight;
%hold on;
cFigure
Y=mesh_struct_IVD2(1).V(:,2); YE=mean(Y(mesh_struct_IVD2(1).E),2);
L=YE<mean(Y);
[Fs,Cs]=element2patch(mesh_struct_IVD2(1).E(L,:),mesh_struct_IVD2(1).C(L));
%subplot(1,2,2);
title('Cut view of the solid tetrahedral mesh model - IVD','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
patch('Faces',Fs,'Vertices',mesh_struct_IVD2(1).V,'FaceColor','flat','CData',Cs,'lineWidth',0.2,'edgeColor',edgeColor);
axis tight;  axis equal;  %grid on;
%colormap(autumn);
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;
%% Mesh vertebrae

fprintf('\n-----------------Solid tetrahedral meshing--------------------\n\n');
%nbodies,stl
stl = 0;
[Lmes] = FU_whole_model_gui(5,stl);


% plot
fontSize =15;

cFigure;
%subplot(1,2,1);

 title('Solid tetrahedral meshing','FontSize',fontSize);
 xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
 patch('Faces',Lmes(1).FT,'Vertices',Lmes(1).VT,'FaceColor','flat','CData',Lmes(1).C,'lineWidth',0.2,'edgeColor',edgeColor);
axis tight;  axis equal;  %grid on; 
camlight headlight;
hold on;
cFigure
Y=Lmes(1).VT(:,1); YE=mean(Y(Lmes(1).E),2);
L=YE<mean(Y);
[Fs,Cs]=element2patch(Lmes(1).E(L,:),Lmes(1).C(L));
%subplot(1,2,2);
title('Cut view of the solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
patch('Faces',Fs,'Vertices',Lmes(1).VT,'FaceColor','flat','CData',Cs,'lineWidth',0.2,'edgeColor',edgeColor);
axis tight;  axis equal;  %grid on;
colormap(autumn);
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% orientation post mesh
[CM, CM_IVD,mesh_struct_IVD2,Lmes] = orientation_gui(Lmes,mesh_struct_IVD2,nbodies, alpha, IVD, hL);

%-----------------------------   plot bodies oriented -------------
%figure
%patch('Faces',Lmes(5).FT,'Vertices',Lmes(5).VT,'FaceColor','r','faceAlpha',0.8,'edgeColor','k','lineWidth',0.5);%'flat','CData',C,'lineWidth',edgeWidth,'edgeColor',edgeColor);

cFigure;
title('3D orientation','FontSize',fontSize);
for(j=1:4)
patch('Faces',mesh_struct_IVD2(j).FE,'Vertices',mesh_struct_IVD2(j).V,'FaceColor','b','faceAlpha',0.8,'edgeColor','k','lineWidth',0.5);%'flat','CData',C,'lineWidth',edgeWidth,'edgeColor',edgeColor);
hold on;
end
for(k=1:5)    
patch('Faces',Lmes(k).FT,'Vertices',Lmes(k).VT,'FaceColor','r','faceAlpha',0.8,'edgeColor','k','lineWidth',0.5);%'flat','CData',C,'lineWidth',edgeWidth,'edgeColor',edgeColor);
hold on;
end
axis tight;  axis equal;  grid on;
colormap(autumn);
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;
%--------------------------------------------------------------


%% write & run FEBio
nbodies = 2;    
fprintf('----------------- FE pre-processing: FEBio --------------------\n\n');
FEB_struct = febio_preprocessing_gui2(febName,nbodies, Lmes,mesh_struct_IVD2);


%% 
%
% 
% 
% _*LMG*_ 
% <https://celavecchia.github.io/LMG/>
% 
% _Carolina Eleonora Lavecchia_, <lavecchia.carolina@gmail.com>