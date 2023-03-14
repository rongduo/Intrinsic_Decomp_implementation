%get the normal map from depth map
function  N=getNormal(Z)
	V=getVectors(size(Z,1),size(Z,2));
	V=V.*Z(:,:,[1 1 1]);
	[nx ny nz]=surfnorm(V(:,:,1),V(:,:,2),V(:,:,3)); %[Nx,Ny,Nz] = surfnorm(...) returns the components of the three-dimensional surface normals for the surface.
	N=cat(3,nx,ny,nz); 
end
