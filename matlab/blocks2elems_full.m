function [xelems,yelems,zelems]=blocks2elems_full(xyblocks,Nx,Col,Row,Lay,ne,z,Height,P)

% Landon Brockmeyer
% 6-27-14

%load save3
%load ess_xyblocks
%Lay=6*Lay;
%Height=HD;
%P=PD;


% Take information in xyblocks and store it in a more convinient array
% Also copy the center pin information to the identical inside pins.

% xyblocks(rad pts, circ pts, elms, 2, ht pts)


% Make it this order because Matlab reads column major order
% xelems(GLL R,GLL C,GLL L,Row,Column,blocks,Layer)

np=3*ne*ne-3*ne+1; % total number of pins
nu=6*ne-5; % number of unique pins
ni=np-nu; % number of inside/non-unique pins
block=size(xyblocks,3); % Only the unique elements 

block_all=block+6*sum(1:ne-2)*6; % Inculding the non-unique inside elements

%xyzelems=zeros(3,Nx+1,Nx+1,Nx+1,Row,Col,elem_all,Lay);
xelems=zeros(Nx+1,Nx+1,Nx+1,Row,Col,block_all,Lay);
yelems=zeros(Nx+1,Nx+1,Nx+1,Row,Col,block_all,Lay);
zelems=zeros(Nx+1,Nx+1,Nx+1,Row,Col,block_all,Lay);
%xyzelems=zeros(Lay,elem_all,Col,Row,Nx+1,Nx+1,Nx+1,3);


% Array that given a Layer and gl, finds corresponding h
L2h=zeros(Lay,Nx+1);
for L=1:Lay
   for gl=1:Nx+1
      L2h(L,gl)=Nx*L-(Nx-gl);
   end
end
% Array that given a Column and gc, finds corresponding c
C2c=zeros(Col,Nx+1);
for C=1:Col
   for gc=1:Nx+1
      C2c(C,gc)=Nx*C-(Nx-gc);
   end
end
% Array that given a Row and gr, finds corresponding r
R2r=zeros(Row,Nx+1);
for R=1:Row
   for gr=1:Nx+1
      R2r(R,gr)=Nx*R-(Nx-gr);
   end
end


% Reorder pins, incrementing left to right then bottom up
% Leave space for the non-unique inside pins

[cpin,epin]=make_arrays(ne);

esspin=find(cpin>0); % position of ess_pins



new_block=zeros(block,1);
trackp=1;
trackv=1;
hole=zeros((block_all-block)/6,2);
hole_n=1;

for i=1:np
   if any(esspin==i)
      new_block(trackp:trackp+4)=trackv:trackv+4;
      trackp=trackp+5;trackv=trackv+5;
      if i==(np+1)/2
         new_block(trackp)=trackv;
         trackp=trackp+1;trackv=trackv+1;
      end
   else
      hole(hole_n,:)=[trackv,trackv+5];
      hole_n=hole_n+1;
      trackv=trackv+6;
   end
end
for i=(trackv+1):block_all
   new_block(trackp)=trackv;
   trackp=trackp+1;trackv=trackv+1;
end
new_block(block)=block_all;

   
  

% Store all of the unique elements in xyzelems

for gr=1:Nx+1
 for gc=1:Nx+1
  for gl=1:Nx+1
   for R=1:Row
    for C=1:Col
     for L=1:Lay
      %xyzelems(1,gr,gc,gl,R,C,new_block,L)=xyblocks(R2r(R,gr),C2c(C,gc),:,1,L2h(L,gl));
      %xyzelems(2,gr,gc,gl,R,C,new_block,L)=xyblocks(R2r(R,gr),C2c(C,gc),:,2,L2h(L,gl));
      %xyzelems(3,gr,gc,gl,R,C,new_block,L)=Height/Lay*(L+z(gl)-1);        
        
      xelems(gr,gc,gl,R,C,new_block,L)=xyblocks(R2r(R,gr),C2c(C,gc),:,1,L2h(L,gl));
      yelems(gr,gc,gl,R,C,new_block,L)=xyblocks(R2r(R,gr),C2c(C,gc),:,2,L2h(L,gl));
      zelems(gr,gc,gl,R,C,new_block,L)=Height/Lay*(L+z(gl)-1);
     end
    end
   end
  end
 end
end

% Copy the center pin information to the identical inside pins

% Find the theta and radius of each hole position xyzelems
xy=zeros(ni,2);
trackxy=1;
for i=1:6
 for j=1:ne-2
  for k=1:j
   %th=(k-1)/j*pi/3+(i-1)*pi/3; %WRONG

   dx=(j-(k-1)/2)*P;
   dy=(k-1)*sqrt(P^2-(P/2)^2);
   th=atan2(dy,dx)+(i-1)*pi/3;

   r=sqrt(dx^2+dy^2);
   xy(trackxy,1:2)=[r*cos(th),r*sin(th)];
   trackxy=trackxy+1;
  end
 end
end



   
% Sort the copied pins by most bottom, then most left

for i=2:ni
   for j=i:-1:2
      if xy(j,2)-xy(j-1,2)<-1e-6
         store=xy(j,:);
         xy(j,:)=xy(j-1,:);
         xy(j-1,:)=store;
      elseif abs(xy(j,2)-xy(j-1,2))<1e-6 && xy(j,1)<xy(j-1,1)
         store=xy(j,:);
         xy(j,:)=xy(j-1,:);
         xy(j-1,:)=store;
      end
   end
end

         


ce=5*(nu-1)/2+6*ni/2; % starting center pin element   

for gr=1:Nx+1
 for gc=1:Nx+1
  for gl=1:Nx+1
   for R=1:Row
    for C=1:Col
     for i=1:ni
      for L=1:Lay
       %xyzelems(1,gr,gc,gl,R,C,hole(i,1):hole(i,2),L)=(1,gr,gc,gl,R,C,ce+(1:6),L)+xy(i,1);
       %xyzelems(2,gr,gc,gl,R,C,hole(i,1):hole(i,2),L)=xyzelems(2,gr,gc,gl,R,C,ce+(1:6),L)+xy(i,2);
       %xyzelems(3,gr,gc,gl,R,C,hole(i,1):hole(i,2),L)=xyzelems(3,gr,gc,gl,R,C,ce+(1:6),L);  
         
       xelems(gr,gc,gl,R,C,hole(i,1):hole(i,2),L)=xelems(gr,gc,gl,R,C,ce+(1:6),L)+xy(i,1);
       yelems(gr,gc,gl,R,C,hole(i,1):hole(i,2),L)=yelems(gr,gc,gl,R,C,ce+(1:6),L)+xy(i,2);
       zelems(gr,gc,gl,R,C,hole(i,1):hole(i,2),L)=zelems(gr,gc,gl,R,C,ce+(1:6),L);            
      end
     end
    end
   end
  end
 end
end



%save ess_xyblocks pin_blcko pin_block type_block xelems yelems zelems ne Dw Df;










