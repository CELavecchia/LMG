function [model, facets,nodes] = obtain_geometry_lam( L_body)

warning 'off'

shp = alphaShape(L_body(:,1),L_body(:,2),L_body(:,3),9);
%plot(shp)
[facets,nodes] = boundaryFacets(shp);
[tri,V] = alphaTriangulation(shp);
%tri = triangulateFaces(shp.tri);
%trimesh(tri,V(:,1),V(:,2),V(:,3));


[N]=numConnect(facets,nodes);
logicThree=N==3;
%remove3connected points
[Ft,Vt,~,L]=triSurfRemoveThreeConnect(facets,nodes,[]);
C=double(L);
%surface smoothening
%{
cPar.n=30;      %check the best value
cPar.Method='HC';
[Vt]=patchSmooth(Ft,Vt,[],cPar);
%}

%get the model
model = createpde();

geometryFromMesh(model,Vt',Ft');

generateMesh(model,'GeometricOrder','linear','Hmax',1.5,'Hmin',1.3);
%pdeplot3D(model), hold on;

%figure;
%pdegplot(model,'FaceAlpha',1), hold on;

%generateMesh(model,'Hmax',1);
%pdeplot3D(model), hold on;

end