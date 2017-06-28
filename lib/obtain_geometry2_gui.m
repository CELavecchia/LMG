function [pdem, facets, nodes] = obtain_geometry2(L_body)%
w = warning ('off','all');
rmpath('./output')
warning(w)
%L = [L_ped_dx;L_body];%lam_new];
L = [L_body];
%splot3(L_ped_dx(:,1),L_ped_dx(:,2),L_ped_dx(:,3),'.r');
shp = alphaShape(L(:,1),L(:,2),L(:,3),2);
%plot(shp);
[facets,nodes] = boundaryFacets(shp);

[tri,V] = alphaTriangulation(shp);
%tri = triangulateFaces(shp.tri);
%trimesh(tri,V(:,1),V(:,2),V(:,3)),hold on;

%{
[N]=numConnect(facets,nodes);
logicThree=N==3;
%remove3connected points
[Ft,Vt,~,L]=triSurfRemoveThreeConnect(facets,nodes,[]);
C=double(L);
%surface smoothening
cPar.n=25;      %check the best value
cPar.Method='HC';
[Vt]=patchSmooth(Ft,Vt,[],cPar);

[v]=tetVolMeanEst(Ft,Vt); 
%}
model = createpde();
%'Hmax',3,'Jiggle','on'
geometryFromMesh(model,nodes',facets');
%generateMesh(model,'GeometricOrder','linear','Hmax',7,'Hmin',0.1);
msh = generateMesh(model,'GeometricOrder','linear','Hmax',1.5,'Hmin',1.0);%it was 7 %'Hmax',7
nodes = msh.Nodes;
elems = msh.Elements;
pdem = createpde;
geometryFromMesh(pdem,nodes,elems);


%pdeplot3D(model), hold on;
%pdegplot(model,'FaceLabels','on','FaceAlpha',1), hold on;

%pdeplot3D(model), hold on;

end