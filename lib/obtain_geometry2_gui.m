function [pdem, facets, nodes] = obtain_geometry2(L_body)%
warning 'off'

%L = [L_ped_dx;L_body];%lam_new];
L = [L_body];
%splot3(L_ped_dx(:,1),L_ped_dx(:,2),L_ped_dx(:,3),'.r');
shp = alphaShape(L(:,1),L(:,2),L(:,3),2);
%plot(shp);
[facets,nodes] = boundaryFacets(shp);

[tri,V] = alphaTriangulation(shp);

model = createpde();
%'Hmax',3,'Jiggle','on'
geometryFromMesh(model,nodes',facets');
%generateMesh(model,'GeometricOrder','linear','Hmax',7,'Hmin',0.1);
msh = generateMesh(model,'GeometricOrder','linear','Hmax',1.5,'Hmin',1.3);%it was 7 %'Hmax',7
nodes = msh.Nodes;
elems = msh.Elements;
pdem = createpde;
geometryFromMesh(pdem,nodes,elems);



end