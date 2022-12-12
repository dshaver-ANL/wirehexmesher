function [xyblocks]=interplevels(xyblocks,m,Nx,Lay,z)

% Landon Brockmeyer
% 6-27-14

% Givin the layers that make up the mesh, interpolates between the
% layers to give the GLL points between layers. 

% xyblocks(circ pts/elem, rad pts/elem, element, coordinates, GLL ht pts)

m1=size(xyblocks,1);
m2=size(xyblocks,2);
ne=size(xyblocks,3);
nH=size(xyblocks,5);


distH=zeros(Lay*Nx+1,1);
n=0;
for i=1:Lay
   for j=1:Nx
      n=n+1;
      distH(n)=(i-1)/(Lay) + z(j)/(Lay);
   end
end
distH(n+1)=1;


for e=1:ne
   for i=1:m1
      xx=spline(1:m,xyblocks(i,1:m2,e,1,1:m),1+(m-1)*distH);
      yy=spline(1:m,xyblocks(i,1:m2,e,2,1:m),1+(m-1)*distH);
      xyblocks(i,1:m2,e,1,:)=xx;
      xyblocks(i,1:m2,e,2,:)=yy;
   end
end


         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         