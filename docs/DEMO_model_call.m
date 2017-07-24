clear; close all; clc;


%% Plotting settings
fontSize=15;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;

%% Control parameters

nbodies = 5;

% file name
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
outputpathName=fullfile(defaultFolder,'data','output','feb');
febName=fullfile(outputpathName,'model_febio');

%%

%function to call the gui and start the model
fprintf('--------------------LMG: LUMBAR MODEL GENERATOR------------------\n\n\n')

%[model, var] =GUI_test1
%{
get:
- type of model (a.average, b.subject-specific)
- num of bodies (you can use the whole spine or only a functional unit.. it
has to be improved)
- output (.stl,.feb,.inp)
- alpha (lumbar curvature angle, if zero then the bodies are arranged in vertical)

%}

%% Build default average geometrical model 
% Provide explanation here

fprintf('-----------------Build the geometrical model--------------------\n');
%  [dimensions] = parameter_average(model.average.age,model.average.height
%  ,model.average.sex); % from output gui

[dimensions] = averageModelInput;

alpha = 43.49;

% fitting and parameterize the model
mesh_struct_IVD2 = model_elaboration_gui(dimensions);

%%
% Visualizing meshes

Y=mesh_struct_IVD2(1).V(:,2); YE=mean(Y(mesh_struct_IVD2(1).E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(mesh_struct_IVD2(1).E(L,:),mesh_struct_IVD2(1).C(L));

cFigure;
subplot(1,2,1);
hold on;
title('Solid hexahedral meshing IVD','FontSize',fontSize); 
gpatch(mesh_struct_IVD2(1).Fb,mesh_struct_IVD2(1).V,mesh_struct_IVD2(1).faceBoundaryMarker);
axisGeom(gca,fontSize); 
camlight headlight;
colormap(gjet);

subplot(1,2,2);
hold on;
title('Cut view of the solid hexahedral mesh model - IVD','FontSize',fontSize);
gpatch(Fs,mesh_struct_IVD2(1).V,Cs);
axisGeom(gca,fontSize); 
camlight headlight;
colormap(gjet);
drawnow;

%% Mesh vertebrae

fprintf('\n-----------------Solid tetrahedral meshing--------------------\n\n');
%nbodies,stl
stl = 0;
[Lmes] = FU_whole_model_gui(5,stl);

%%
% Visualizing 

Y=Lmes(1).VT(:,1); YE=mean(Y(Lmes(1).E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(Lmes(1).E(L,:),Lmes(1).C(L));

cFigure;
%subplot(1,2,1);
hold on; 
title('Solid tetrahedral meshing','FontSize',fontSize);
patch('Faces',Lmes(1).FT,'Vertices',Lmes(1).VT,'FaceColor','flat','CData',Lmes(1).C,'lineWidth',0.2,'edgeColor',edgeColor);
axis tight;  axis equal;  %grid on;
colormap(autumn);
axisGeom(gca,fontSize); 
camlight headlight;
drawnow;

cFigure;
hold on; 

%subplot(1,2,2);
title('Cut view of the solid tetrahedral mesh model','FontSize',fontSize);
patch('Faces',Fs,'Vertices',Lmes(1).VT,'FaceColor','flat','CData',Cs,'lineWidth',0.2,'edgeColor',edgeColor);
axisGeom(gca,fontSize); 
colormap(autumn);
camlight headlight;
drawnow;

%% orientation post mesh

[CM, CM_IVD,mesh_struct_IVD2,Lmes] = orientation_gui(Lmes,mesh_struct_IVD2, alpha, dimensions.IVD, dimensions.hL);

%% 
% Visualizing ...

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

axisGeom(gca,fontSize);
colormap(autumn);
camlight headlight;
drawnow;

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