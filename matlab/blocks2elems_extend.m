function [xelems,yelems,zelems]=blocks2elems_extend(xyblocks,Nx,Col,Row,Lay,ne,z,Height,P)

% Landon Brockmeyer
% 6-27-14

%load save3
%load ess_xyblocks
%Lay=6*Lay;
%Height=HD;
%P=PD;

nex=ne;
ne=3;


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

%

%xyzelems=zeros(3,Nx+1,Nx+1,Nx+1,Row,Col,elem_all,Lay);
xelemst=zeros(Nx+1,Nx+1,Nx+1,Row,Col,block_all,Lay);
yelemst=zeros(Nx+1,Nx+1,Nx+1,Row,Col,block_all,Lay);
zelemst=zeros(Nx+1,Nx+1,Nx+1,Row,Col,block_all,Lay);
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
new_block;
   
  

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
        
      xelemst(gr,gc,gl,R,C,new_block,L)=xyblocks(R2r(R,gr),C2c(C,gc),:,1,L2h(L,gl));
      yelemst(gr,gc,gl,R,C,new_block,L)=xyblocks(R2r(R,gr),C2c(C,gc),:,2,L2h(L,gl));
      zelemst(gr,gc,gl,R,C,new_block,L)=Height/Lay*(L+z(gl)-1);
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
         
       xelemst(gr,gc,gl,R,C,hole(i,1):hole(i,2),L)=xelemst(gr,gc,gl,R,C,ce+(1:6),L)+xy(i,1);
       yelemst(gr,gc,gl,R,C,hole(i,1):hole(i,2),L)=yelemst(gr,gc,gl,R,C,ce+(1:6),L)+xy(i,2);
       zelemst(gr,gc,gl,R,C,hole(i,1):hole(i,2),L)=zelemst(gr,gc,gl,R,C,ce+(1:6),L);   
       
      end
     end
    end
   end
  end
 end
end
%}





if nex==3
    xelems=xelemst;
    yelems=yelemst;
    zelems=zelemst;
    
    
else
    
%Added 7-8-15 
%Run ineri with ne=3 to cut down time, then turn the ne=3
%bundle into ne=x bundle

% Modify the 19pin case to 37, 61, etc
% There are 13 types of pins, make 19 long array of pin type corresponding
% to position, then do same for the desired, larger size. Pins are numbered
% left to right bottom to top.

npx=3*nex*nex-3*nex+1; % total number of pins nex
nb=6*(np-(ne-1))+(ne*6);
nbx=6*(npx-(nex-1))+(nex*6);
b2t19=zeros(nb,1);
b2tx=zeros(nbx,1);
bnum=zeros(nbx,1);
b2p19=zeros(nb,1);
b2px=zeros(nbx,1);

xelems=zeros(Nx+1,Nx+1,Nx+1,Row,Col,nbx,Lay);
yelems=zeros(Nx+1,Nx+1,Nx+1,Row,Col,nbx,Lay);
zelems=zeros(Nx+1,Nx+1,Nx+1,Row,Col,nbx,Lay);


outsidepinorder=[1 2 3 7 12 16 19 18 17 13 8 4]; % can make
orderoutside=[2 7 16 18 13 4 1 1 3 3 12 12 19 19 17 17 8 8]; %ospo(op19)
outpin19=[2 4 6 8 10 12 1 1 3 3 5 5 7 7 9 9 11 11]';% make easily: order presented around outside bo
outtype19=[ 7 9 11 12 10 8 1 1 2 2 4 4 6 6 5 5 3 3]';% pin types around outside pt19(orderoutside)


% pin type order around outside
pt19=[ 1 7 2 8 13 13 9 3 13 13 13 4 10 13 13 11 5 12 6]';

ptnew=zeros(npx,1);
ind=1;

ptnew(ind)=1; ind=ind+1;
ptnew(ind:ind+nex-2)=7; ind=ind+nex-2;
ptnew(ind)=2; ind=ind+1;

for i=1:nex-2
    ptnew(ind)=8; ind=ind+1;
    ptnew(ind:ind+i+nex-2)=13; ind=ind+i+nex-2;
    ptnew(ind)=9; ind=ind+1;
end

ptnew(ind)=3; ind =ind+1;
ptnew(ind:ind+2*nex-3)=13; ind=ind+2*nex-3;
ptnew(ind)=4; ind=ind+1;

for i=nex-2:-1:1
    ptnew(ind)=10; ind=ind+1;
    ptnew(ind:ind+i+nex-2)=13; ind=ind+i+nex-2;
    ptnew(ind)=11; ind=ind+1;
end

ptnew(ind)=5; ind=ind+1;
ptnew(ind:ind+nex-2)=12; ind=ind+nex-2;
ptnew(ind)=6; 

%boarder blocks outside order
boo19=[2 4 6 8 10 12 1 1 3 3 5 5 7 7 9 9 11 11]';
%{
boox=zeros(6*nex,1);
for i=1:6
    boox((nex-2)*(i-1)+1:(nex-2)*i)=boo19(i);
end
boox((nex-2)*6:end)=boo19((ne-2)*6:end);
%}
boox=zeros(6*(nex-1),1);
ind=1;
for i=1:6
    boox(ind:ind+(nex-3))=(nex-1)*i-(nex-3):(nex-1)*i;
    ind=ind+(nex-2);
end

for i=1:6
    boox(ind:ind+1)=(nex-1)*i-(nex-2);
    ind=ind+2;
end

%boarder pins order
bpo19=[1 2 3 7 12 16 19 18 17 13 8 4]';
bpox = zeros(6*(ne-1),1);
bpox(1:nex)=1:nex;
ind=nex+1;
for i=2:size(ptnew)
    if ptnew(i-1)==13 && ptnew(i)~=13
        bpox(ind)=i;
        ind=ind+1;
    end 
end
bpox(ind:ind+nex-1)=size(ptnew):-1:size(ptnew)-nex+1;
ind=ind+nex;
for i=size(ptnew)-1:-1:1
    if ptnew(i+1)==13 && ptnew(i)~=13
        bpox(ind)=i;
        ind=ind+1;
    end 
end


% boarder block insert order by pin
%biop19=[2 7 16 18 13 4 1 1 3 3 12 12 19 19 17 17 8 8];
biop19=bpo19(boo19);
biopx=bpox(boox);


%pt19[ 1 7 2 8 13 13 9 3 13 13 13 4 10 13 13 11 5 12 6]';

% boarder block insert order by type
%biot19=[ 7 9 11 12 10 8 1 1 2 2 4 4 6 6 5 5 3 3]'
biot19=pt19(biop19);
biotx=ptnew(biopx);


% given block find pin

ind=1;
for i=1:size(pt19)
    if pt19(i)~=13
      b2t19(ind:ind+4)=pt19(i); 
      b2p19(ind:ind+4)=i;
      ind=ind+5;
      
    else
      b2t19(ind:ind+5)=pt19(i); 
      b2p19(ind:ind+5)=i;
      ind=ind+6;
    end
end
%{
b2t19
biot19

size(b2t19)
size(b2t19==0)
size(b2t19(ind:end))
pause
%}
b2t19(ind:end)=biot19;
b2p19(ind:end)=biop19;
%given block find pin

ind=1;
for i=1:size(ptnew)
    if ptnew(i)~=13
      b2tx(ind:ind+4)=ptnew(i);
      b2px(ind:ind+4)=i; 
      bnum(ind:ind+4)=1:5;
      ind=ind+5;
    else
      b2tx(ind:ind+5)=ptnew(i);
      b2px(ind:ind+5)=i; 
      bnum(ind:ind+5)=1:6;
      ind=ind+6;
    end
end

b2tx(ind:end)=biotx;
b2px(ind:end)=biopx;
bnum(ind:end)=6;
array=0:2:10;
bnum(end-array)=7;

% Find pin center corresponding to pin number for 19 and 61
    %19
    pne=zeros(np,2);
    trackxy=1;
    for i=1:6
     for j=1:ne-1
      for k=1:j
       %th=(k-1)/j*pi/3+(i-1)*pi/3; %WRONG

       dx=(j-(k-1)/2)*P;
       dy=(k-1)*sqrt(P^2-(P/2)^2);
       th=atan2(dy,dx)+(i-1)*pi/3;

       r=sqrt(dx^2+dy^2);
       pne(trackxy,1:2)=[r*cos(th),r*sin(th)];
       trackxy=trackxy+1;
      end
     end
    end
    pne(end,:)=[0 0];
    % Sort the copied pins by most bottom, then most left
    for i=2:np
       for j=i:-1:2
          if pne(j,2)-pne(j-1,2)<-1e-6
             store=pne(j,:);
             pne(j,:)=pne(j-1,:);
             pne(j-1,:)=store;
          elseif abs(pne(j,2)-pne(j-1,2))<1e-6 && pne(j,1)<pne(j-1,1)
             store=pne(j,:);
             pne(j,:)=pne(j-1,:);
             pne(j-1,:)=store;
          end
       end
    end

    % nex (61)

    pnex=zeros(npx,2);
    trackxy=1;
    for i=1:6
     for j=1:nex-1
      for k=1:j
       %th=(k-1)/j*pi/3+(i-1)*pi/3; %WRONG

       dx=(j-(k-1)/2)*P;
       dy=(k-1)*sqrt(P^2-(P/2)^2);
       th=atan2(dy,dx)+(i-1)*pi/3;

       r=sqrt(dx^2+dy^2);
       pnex(trackxy,1:2)=[r*cos(th),r*sin(th)];
       trackxy=trackxy+1;
      end
     end
    end
    pnex(end,:)=[0 0];
    % Sort the copied pins by most bottom, then most left
    for i=2:npx
       for j=i:-1:2
          if pnex(j,2)-pnex(j-1,2)<-1e-6
             store=pnex(j,:);
             pnex(j,:)=pnex(j-1,:);
             pnex(j-1,:)=store;
          elseif abs(pnex(j,2)-pnex(j-1,2))<1e-6 && pnex(j,1)<pnex(j-1,1)
             store=pnex(j,:);
             pnex(j,:)=pnex(j-1,:);
             pnex(j-1,:)=store;
          end
       end
    end

    %xelems=zeros(Nx+1,Nx+1,Nx+1,Row,Col,nbx,Lay);
    clf
hold on
for i=1:nbx
    % correct pin
    pin=find(b2t19==b2tx(i));
    block=bnum(i);
    j=pin(block);
    xelems(:,:,:,:,:,i,:)=xelemst(:,:,:,:,:,j,:)-pne(b2p19(j),1)+pnex(b2px(i),1);
    yelems(:,:,:,:,:,i,:)=yelemst(:,:,:,:,:,j,:)-pne(b2p19(j),2)+pnex(b2px(i),2);
    zelems(:,:,:,:,:,i,:)=zelemst(:,:,:,:,:,j,:);
    
    %{
    % NOW PUT IN THE  DIFFERENCE IN POSITION
    for gr=1:2:Nx+1
     for gc=1:2:Nx+1
      %for gl=1:Nx+1
       for R=1:2:Row
        for C=1:2:Col
    
            plot(xelems(gr,gc,1,R,C,i,1),yelems(gr,gc,1,R,C,i,1),'b.')
        end
       end
      %end
     end
    end
  %}
end
%pause

end

%save 61pin_coarse5 pin_blcko pin_block type_block xelems yelems zelems ne Dw Df;





