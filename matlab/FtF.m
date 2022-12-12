           %       load -mat 61pin_v5.mat;
 
                  nlay=size(xelems,7); 
                  ncolumn=size(xelems,5);
                  %nrow=size(xelems,6);
                  nrow=size(xelems,4);

                %  xelems_lay=zeros(3,3,3,12,6*ne,nlay);
                %  yelems_lay=zeros(3,3,3,12,6*ne,nlay);
                %  zelems_lay=zeros(3,3,3,12,6*ne,nlay);
                  
                  xelems_lay=zeros(3,3,3,ncolumn,6*ne,nlay);
                  yelems_lay=zeros(3,3,3,ncolumn,6*ne,nlay);
                  zelems_lay=zeros(3,3,3,ncolumn,6*ne,nlay);
                  
                  block1=size(xelems,6)-6*ne;
                  
                  Old_FtF=max(max(max(max(max(max(max(yelems)))))))*2.0;

% or in meters, .1528852416 and .1765366506 m
% while in specs they are 0.15382 and 0.177616 m

                  New_FtF=Old_FtF*FtF_rescale;
                  
                  ff0=New_FtF/Old_FtF;
                  ff1=1.0+(ff0-1.0)*0.5;
                  
                  for lay=1:nlay;
                  for gl=1:3;   
                  for i=1:6*ne; 
                  for j=1:ncolumn;   
                  for k=1:3;     
                  xelems_lay(1,k,gl,j,i,lay)=xelems(3,k,gl,nrow,j,i+block1,lay);
                  yelems_lay(1,k,gl,j,i,lay)=yelems(3,k,gl,nrow,j,i+block1,lay); 
                  zelems_lay(1,k,gl,j,i,lay)=zelems(3,k,gl,nrow,j,i+block1,lay);  
                  xelems_lay(2,k,gl,j,i,lay)=xelems(3,k,gl,nrow,j,i+block1,lay)*ff1;
                  yelems_lay(2,k,gl,j,i,lay)=yelems(3,k,gl,nrow,j,i+block1,lay)*ff1; 
                  zelems_lay(2,k,gl,j,i,lay)=zelems(3,k,gl,nrow,j,i+block1,lay);       
                  xelems_lay(3,k,gl,j,i,lay)=xelems(3,k,gl,nrow,j,i+block1,lay)*ff0;
                  yelems_lay(3,k,gl,j,i,lay)=yelems(3,k,gl,nrow,j,i+block1,lay)*ff0; 
                  zelems_lay(3,k,gl,j,i,lay)=zelems(3,k,gl,nrow,j,i+block1,lay); 
                  end;
                  end;
                  end;
                  end;
                  end;
                  
%                   for i=1:6*ne; 
%                   for j=1:12;   
%                   for k=1:3;
%                       xx=xelems_lay(3,k,1,j,i,1);
%                       yy=yelems_lay(3,k,1,j,i,1);
%                       plot(xx,yy,'o'); hold on;
%                       xx=xelems_lay(2,k,1,j,i,1);
%                       yy=yelems_lay(2,k,1,j,i,1);
%                       plot(xx,yy,'v'); hold on;                   
%                       xx=xelems_lay(1,k,1,j,i,1);
%                       yy=yelems_lay(1,k,1,j,i,1);
%                       plot(xx,yy,'+'); hold on;              
%                   end;
%                   end;
%                   end;

fprintf('begin nekify \n')
tic
[xnek,ynek,znek,BCs,BCcon]=nekify2(xelems,yelems,zelems,xelems_lay,yelems_lay,zelems_lay,ne);
r_data(1)=size(xelems,7);
r_data(2)=size(xelems,6);
r_data(3)=size(xelems,5)*size(xelems,4);
r_data(4)=6*ne*size(xelems,5);
r_data(5)=rperiodic;
save nekinfo r_data xnek ynek znek BCs BCcon
fprintf('nekify complete \n')
toc

fprintf('file contains %i gll points in %i elements \n', numel(xnek)/3, size(xnek,1))
fprintf('FtF layering complete, now printnek \n')
printnek;                  
fprintf('Printing Done! \n')
                  
