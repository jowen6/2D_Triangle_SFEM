function [node,elem,surfacedata] = GenerateMesh(Geo,NumRefinements)

%Generate Mesh
%==========================================================================
switch Geo
    case 1  %Sphere
        surfacedata = spheresurface();
        
    case 2  %Torus
        surfacedata = torussurface_J();
          
    case 3  %ellipsoid
        surfacedata = ellipsoidsurface_J();
        
    case 4  %Quartics  
        surfacedata = quarticssurface();
        
    case 5  %OrthoCircle
        surfacedata = orthocirclesurface();
    
    case 6  %Heart
        surfacedata = heartsurface_J();
        
end

[node,elem] = surfacedata.initmesh();
node = surfacedata.project(node);
%[node,elem] = optsurfacemesh(node,elem,surfacedata);

if NumRefinements > 1
    for i = 1:NumRefinements-1
        [node,elem] = smeshuniformrefine(node,elem);
        node = surfacedata.project(node);
        %[node,elem] = optsurfacemesh(node,elem,surfacedata);
    end
end
%==========================================================================