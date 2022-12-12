% Mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D    = 8.00;            %  pin diameter  (mm)
P    = 9.04;            %  pin pitch
Dw   = 1.03;            %  wire diameter (mm)
Df   = 0.4;             %  fillet diameter (mm) (making this too small can cause bad elements)
H    = 100.0;           %  wire pitch (mm)
T    = 0.05*Dw          % trimmed off of wire (mm) (cuts off the tip of the wire)
S    = 0.025*Dw;        % Wire submerged (mm) (how sunken in the wire is into the pin)
Adjust = 1;             % if 1, Adjust flattness of wire when away from pins (if trimmed is on, only trim when wire passes pin)
iFtF   = 0;             % if 1, add layer next to outer can
G  = 0.0525;            % gap between wire and wall (mm)
FtF_rescale = 1.0086;   % rescaling of outer FtF, the difference is given by an additional mesh layer - MUST BE BIGGER THAN 1
ne=2;                   % Number of edge pins. e.g., ne=3 for 19 pin assembly. MINIMUM is 2
Col=12;                 % Number of columns per block (5 or 6 blocks surround each pin)
Row=3;                  % Number of rows per block (2 blocks fit between neighboring pins)
Rowdist=[65 30 5];      % distribution of elements in the row (for creating boundary layer)
Lay=20;                 % Number of layers for 60 degree turn (this should be a reasonable number - above 10, tested typically for ~20)
rperiodic=1.0;          % Less than zero if periodic, Greater than zero if inlet/outlet
ipolar=0.25*(D/H)*360   % Starting polar orientation of wire (in degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PD =  P/D; 
Dw = Dw/D; 
Df = Df/D; 
GD = G/D;
HD =  H/D;
TD=T/D;
SD=S/D;
a=Adjust;

Nx=2;    % Number of GLL integration points - DO NOT CHANGE THIS LINE OR THE FOLLOWING!!
Plane=3;
         
[z,wt]=zwgll(Nx); % z is distribution of gll points from -1 to 1
z=z/2+0.5; % change distribution to 0 to 1       
Rowdist=Rowdist/100*Row; 

tstart=tic;
tic
fprintf('begin wire_mesher \n')
fprintf('begin build_wire_template \n')
build_wire_template_L(PD,Dw,Df,TD,SD)   %Uncomment on first run with geometry
if a==1
    build_wire_template_L(PD,Dw,Df,0,SD)
end
fprintf('build_wire template complete \n')
toc


m=Lay*(Plane-1)+1;   % Number of planes
dth=(pi/3)/(Lay*(Plane-1));
th_start=(pi/180.0)*ipolar; %Radians
mg=6*(Lay)*Nx+1;% Number of verital gll points

%all_blocks = 0;  % If you want blocks only for center pin; NOT VALIDATED (11/19/09)
all_blocks = 1;  % If you want blocks for all 13 pins

tic


[pin_blcko,pin_block,type_block,xyblock]=gen_ess_blocks_L(th_start,0,Nx,Col,Row,Rowdist,z,3,TD,SD,a,GD); hold off; 
%[pin_blcko,pin_block,type_block,xyblock]=gen_ess_blocks_L(pi/6,0,Nx,Col,Row,Rowdist,z,3,TD,SD,a); hold off;
fprintf('layer 1/%i complete \n',m)
toc 

s1=size(xyblock,1); s2=size(xyblock,2); s3=size(xyblock,3); s4=size(xyblock,4); 

xyblocks=zeros(s1,s2,s3,s4,mg); 

kk=1; xyblocks(:,:,:,:,kk)=xyblock;


for k=1:m-1  ; th=k*dth+th_start; kk=kk+1; 
     tic
   %[pin_blcko,pin_block,type_block,xyblocks(:,:,:,:,kk)]=gen_ess_blocks_L(th,k,Nx,Col,Row,Rowdist,z,ne,TD,SD);extend
[pin_blcko,pin_block,type_block,xyblocks(:,:,:,:,kk)]=gen_ess_blocks_L(th,k,Nx,Col,Row,Rowdist,z,3,TD,SD,a,GD);   
   hold off; 
    fprintf('layer %d/%i complete \n', k+1,m)
   toc
 end;
 
save saverotate xyblocks pin_block ne
 
 
% Takes first 60 degrees and copies it to create a full rotation 
fprintf('begin rotate \n')
tic
for i=1:5
   lprev=(i-1)*(m-1)+1; 
   lshare=i*(m-1)+1; 
   lnext=(i+1)*(m-1)+1; 
   %xyblocks(:,:,:,:,lshare:lnext)=rotate(xyblocks(:,:,:,:,lprev:lshare),pin_block,ne);extend
      xyblocks(:,:,:,:,lshare:lnext)=rotate(xyblocks(:,:,:,:,lprev:lshare),pin_block,3);
end
fprintf('rotate complete \n')
toc

save save2 xyblocks m Nx Lay z

% Interpolates between existing planes to find GLL points
fprintf('begin interplevels \n')
tic
[xyblocks]=interplevels(xyblocks,6*(m-1)+1,Nx,6*Lay,z);
fprintf('interplevels complete \n')
toc

save save3 xyblocks Nx Col Row Lay ne z HD PD pin_block pin_blcko
%pause
% Takes all of the points and organizes them by GLL element
fprintf('begin blocks2elems_full \n')
tic
%[xelems,yelems,zelems]=blocks2elems_extend(xyblocks,Nx,Col,Row,6*Lay,ne,z,HD,PD);extend

[xelems,yelems,zelems]=blocks2elems_extend(xyblocks,Nx,Col,Row,6*Lay,ne,z,HD,PD);

fprintf('blocks2elems_extend complete \n')
toc
%pause
save 61pin_v5 pin_blcko pin_block type_block xelems yelems zelems ne Dw Df rperiodic T PD;

if (iFtF==1)
FtF;
else

% Reorganizes information for nek, figures BCs
fprintf('begin nekify \n')
[xnek,ynek,znek,BCs,BCcon]=nekify(xelems,yelems,zelems,ne);
r_data(1)=size(xelems,7);
r_data(2)=size(xelems,6);
r_data(3)=size(xelems,5)*size(xelems,4);
r_data(4)=0.0;
r_data(5)=rperiodic;
save nekinfo r_data xnek ynek znek BCs BCcon
fprintf('nekify complete \n');
fprintf('file contains %i gll points in %i elements \n', numel(xnek)/3, (36*ne-24)*Row*Col*6*Lay);
fprintf('ineri complete \n');
printnek;
fprintf('Printing Done! \n');
%fprintf('Solid mesh being constructed ... \n');
%pinmesh5;
%fprintf('Printing Solid mesh ...\n');
%printsolid;
%fprintf('Printing Done! \n');
end;
exit;
