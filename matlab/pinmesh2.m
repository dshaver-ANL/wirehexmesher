% pinmesh_simple.m

% Landon Brockmeyer
% 8-14-14

clf 
hold on
axis equal

tic
fprintf('begin pinmesh_simple \n')

load -mat 61pin_v5.mat  % pin_blcko pin_block type_block xelems yelems zelems ne Dw Df rperiodic T;
%P=8.40;
%D=6.55;

eps_d2=5*1e-6;
%PD=P/D;
%ne=3;

%=================== PARAMETERS, inner dimater Di, outer Diameter Do ==================
Do=1.0;
Di=5.07/6.55;
Dwo=Dw;%1.7/6.55;
Dwi=0.12;%Dwo-(Do-Di);
Dwi=0.8*Dwo;
%Nx=4;
%======================================================================================

[npin,X,Y,xm,ym,km]=pin_coors(ne,PD);

%plot(X,Y,'mx')
%plot(xm,ym,'bx')

npins=npin;
Nx=size(xelems,1)-1;
Col=size(xelems,5);
Lay=size(xelems,7);
Laydo=Lay;
nblocks=size(xelems,6)-6*ne;
%pause

[z,w]=zwgll(Nx);
z1=z;
z=z/2+.5;

     for i=1:Nx+1; 
      P1   = 1.0;
      P2   = 2*z(i)-1;
      P3   = P2;
     for k=1:(Nx-1) 
         Fk  = k;
         P3  = ((2.*Fk+1.)*(2*z(i)-1)*P2 - Fk*P1)/(Fk+1.);
         P1  = P2;
         P2  = P3;
     end;
    plgl(i)=P3;
    end; 
     
      Fn = Nx;
      d0 = Fn*(Fn+1.0)/4.;
      for i=1:Nx+1;
      for j=1:Nx+1;
         Dt(i,j) = 0.0;
         if (i==j)
             if (i==1)      
             Dt(i,i) = -d0;
             end
             if (i==(Nx+1))
             Dt(i,j) =  d0;   
             end;
         else 
         xlgl=plgl(i);
         ylgl=plgl(j);
         Dt(i,j) = xlgl/(ylgl*(z1(i)-z1(j)));
         end;     
      end;
      end;

pinmap=zeros(nblocks,1);

% Important angles
dRW = (Dwo+Do)/2;
dRF = (Do+Df)/2;
dWF = (Dwo+Df)/2;

aa = acos((dRF^2 + dWF^2 - dRW^2)/(2*dRF*dWF));
ab = acos((dRW^2 + dWF^2 - dRF^2)/(2*dRW*dWF));
ac = acos((dRW^2 + dRF^2 - dWF^2)/(2*dRW*dRF));

pmax=ac; 
pmin=-ac; 
pmin2=2*pi-ac;
wmax=-pi+ab;
wmin=pi-ab;
fnpin=-ac+pi;
fnwir=-ac+pi-aa;
fppin=ac-pi;
fpwir=ac-pi+aa;

bmin=-pi/3;
bmax=pi/3;


% Wire and fillet locations when pin centered at 0 and th=0
Xwm=(Do+Dwo)/2;
Ywm=0;
%plot(Xwm,Ywm,'m*')

Xfnm=(Do+Df)/2*cos(-ac);
Yfnm=(Do+Df)/2*sin(-ac);
%plot(Xfnm,Yfnm,'k*')

Xfpm=(Do+Df)/2*cos(ac);
Yfpm=(Do+Df)/2*sin(ac);
%plot(Xfpm,Yfpm,'k*')
%pause

% Define inside and outside boundaries
fine=200;
% fine=100;
%tfine= ??? *fine;




R=1; % Number of rows
% pmesh gives the x,y,and z coordinates for each gll point
% pmesh(Layers, number of blocks, number of columns, number of rows, GLL C, GLL R, coordinates)
%pmesh=zeros(Lay,nblocks*Col,R,Nx+1,Nx+1,Nx+1,3);

levelem=nblocks*Col*R;
totelem=levelem*Laydo;
clear r_data;
r_data(1)=levelem;
r_data(2)=totelem;

% pmesh lists
%pmeshx(totelem,gll r, gll c,gll l)
pmeshx=zeros(totelem,Nx+1,Nx+1,Nx+1);
pmeshy=pmeshx;
pmeshz=pmeshx;

% pmesh_BCs gives the BC for each element's face (0 = no BC, 1 = heat flux, 2 = periodic, 3 = inside wire)
% pmesh_BCs(totelem, 6 sides)    
% sides: 1=inward, 2=cw, 3=outward, 4=ccw, 5=down, 6=up
pmesh_BCs=zeros(totelem,6);

% The elements don't form a perfect circle inside, so heat flux must be adjusted.  
% pmesh_BCi gives the intensity multiplier for the heat flux for each element.
pmesh_BCi=zeros(totelem,1);

% pmesh_BCc gives the element corresponding to each periodic BC
pmesh_BCc=zeros(totelem,1);

%ne=5;
NRR=(2*(ne-1)+1);
i1=-1;
p=0;
for i0=1:NRR;
    if(i0<(ne+1)) i1=i1+1;     end;   
    if(i0>(ne)) i1=i1-1;  end; 
    
    for j1=1:(ne+i1);        
    p=p+1;
    pp(i0,j1)=p;
    pp1(i0,j1)=5;
    bb(p)=5;
    if ((j1>1)&&(j1<(ne+i1)))
        if ((i0>1)&&(i0<NRR))
        pp1(i0,j1)=6;
        bb(p)=6;
        end;
    end;    
    end;    

end;
%bb

size(bb,2)
track=1;
for i=1:size(bb,2)
    pinmap(track:track+bb(i)-1)=i;
    track=track+bb(i);
end
%pinmap

tic

imp_pins=[1 3 5 12 18 27 31 35 44 50 57 59 61];

count =0
for level=1:Laydo
    level

for p=1:npins   
    for gl=1:2;
         
         %clf 
         %hold on
         %axis equal

         con1=0;
         con2=0;
         c2w=0;
         w2c=0;

         block1=find(pinmap==p,1,'first');
         blockf=find(pinmap==p,1,'last'); 
         
         nbb=bb(p);
         nelp=nbb*Col;   


         if gl==1
         count = count+nelp
         level
         p
         end 

         npts=nbb*Col*2+1; 

         bound_out=zeros(npts,2);
         bound_in=zeros(npts,2);

         for i=1:nbb
            for j=1:Col
               for k=1:Nx+1
                  kk=k+(j-1)*2+Col*2*(i-1);
                  bound_out(kk,1)=xelems(1,k,gl,1,j,i+block1-1,level);
                  bound_out(kk,2)=yelems(1,k,gl,1,j,i+block1-1,level); 
                  
               end
            end
           
         end
         
         type=zeros(npts);
       
         %plot(bound_out(:,1),bound_out(:,2),'mx')
         %pause
                 
               for i=1:(npts)
    
                   if (i<npts)
                    bdx=bound_out(i+1,1)-bound_out(i,1);
                    bdy=bound_out(i+1,2)-bound_out(i,2);
                   end
                   if (i==npts)
                    bdx=bound_out(2,1)-bound_out(i,1);
                    bdy=bound_out(2,2)-bound_out(i,2);
                   end
                   
                  
                   
                 th=(level-1+z(gl))*2*pi/Lay;

                 if th>=pi;th=th-2*pi;end

                 Xwm=X(p)+(Do+Dwo)/2*cos(th);
                 Ywm=Y(p)+(Do+Dwo)/2*sin(th);
         %plot(Xwm,Ywm,'k*')

                 Xfnm=X(p)+(Do+Df)/2*cos(-ac+th);
                 Yfnm=Y(p)+(Do+Df)/2*sin(-ac+th);
         %plot(Xfnm,Yfnm,'k*')

                 Xfpm=X(p)+(Do+Df)/2*cos(ac+th);
                 Yfpm=Y(p)+(Do+Df)/2*sin(ac+th);
         %plot(Xfpm,Yfpm,'k*')


               angle_p=atan2(bound_out(i,2)-Y(p),bound_out(i,1)-X(p));
               radius=sqrt((bound_out(i,1)-X(p))^2+(bound_out(i,2)-Y(p))^2);

               if abs(th-angle_p)>pi
                   if th<angle_p;angle_p=angle_p-2*pi;
                   else angle_p=angle_p+2*pi;
                   end
               end


%%%%%%%%%% NEW MESHER
          
%        figure(1) 
%        hold on
         
         bo_x = bound_out(i,1);
         bo_y = bound_out(i,2);


         %distance from wire center
         %distance from fillet 1
         df1 = ((bo_x-Xfnm).^2+(bo_y-Yfnm).^2).^0.5;
         dw = ((bo_x-Xwm).^2+(bo_y-Ywm).^2).^0.5;
          %distance from fillet 2
         df2 = ((bo_x-Xfpm).^2+(bo_y-Yfpm).^2).^0.5;
     
         % fillet bias angle
         fba = 7*pi/12;
 

         if (dw<=Dwo/2*1.01)
    
           % angle to wire center
           aw=atan2(bound_out(i,2)-Ywm,bound_out(i,1)-Xwm);

           bound_in(i,1) = bound_out(i,1)-((Do-Di)/4.0)*cos(aw);
           bound_in(i,2) = bound_out(i,2)-((Do-Di)/4.0)*sin(aw);
 
%           figure(1)
%           plot(bound_out(i,1),bound_out(i,2),'m*')
%           plot(bound_in(i,1),bound_in(i,2),'mx')
%           pause
        
        else

         % inside boundary is x distance toward center of pin from outside boundary
         bound_in(i,1)=bound_out(i,1)-((Do-Di)/4.0)*cos(angle_p);
         bound_in(i,2)=bound_out(i,2)-((Do-Di)/4.0)*sin(angle_p);
 
%           figure(1)
%           plot(bound_out(i,1),bound_out(i,2),'k*')
%           plot(bound_in(i,1),bound_in(i,2),'kx')
%           pause
 
        end

         if (df1<=7*Df/2)
           
            % Need to tranistion from pin and wire boundary to fillet
            % Introduce a factor based on distance from center of fillet

            factor = (df1/(Df/2)-1)/6;

            %angle to fillet 1

            %bias to a certain angle
           af1=atan2(bound_out(i,2)-Yfnm,bound_out(i,1)-Xfnm);

           aoi = [th+fba th-(2*pi-fba) th-(4*pi-fba) th+2*pi+fba];
           iii = find(abs(aoi-af1)<pi);
 
           if(size(iii,2)~=1)
              af1
              iii
              aoi
              pause
           end

          af1 = (af1+2*aoi(iii))/3;

           bin_x = bound_out(i,1)+((Do-Di)/4.0)*cos(af1);
           bin_y = bound_out(i,2)+((Do-Di)/4.0)*sin(af1);

           % Transition using factor
           bound_in(i,1) = bin_x*(1-factor)+bound_in(i,1)*factor;
           bound_in(i,2) = bin_y*(1-factor)+bound_in(i,2)*factor;
         

%           figure(1)
%           plot(bound_out(i,1),bound_out(i,2),'g*')
%           plot(bound_in(i,1),bound_in(i,2),'gx')
%           pause
        
         elseif (df2<=7*Df/2)
           
            % Need to tranistion from pin and wire boundary to fillet
            % Introduce a factor based on distance from center of fillet

            factor = (df2/(Df/2)-1)/6;


            %bias to a certain angle

           af2=atan2(bound_out(i,2)-Yfpm,bound_out(i,1)-Xfpm);
           aoi = [th-fba th+(2*pi-fba) th+(4*pi-fba) th-2*pi-fba];
           iii = find(abs(aoi-af2)<pi);

           if(size(iii,2)~=1)
              af2
              iii
              aoi
              pause
           end

           af2 = (af2+2*aoi(iii))/3;

           bin_x = bound_out(i,1)+((Do-Di)/4.0)*cos(af2);
           bin_y = bound_out(i,2)+((Do-Di)/4.0)*sin(af2);

           % Transition using factor
           bound_in(i,1) = bin_x*(1-factor)+bound_in(i,1)*factor;
           bound_in(i,2) = bin_y*(1-factor)+bound_in(i,2)*factor;

%{
           bound_in(i,1) = bound_out(i,1)+((Do-Di)/4.0)*cos(af2);
           bound_in(i,2) = bound_out(i,2)+((Do-Di)/4.0)*sin(af2);

           % bias to a certain point
           poix = X(p) + Do/2*cos(th); 
           poiy = Y(p) + Do/2*sin(th);
           bin_x = (bound_in(i,1)+poix)/2;
           bin_y = (bound_in(i,2)+poiy)/2;
%}

           % Transition using factor
%           bound_in(i,1) = bin_x*(1-factor)+bound_in(i,1)*factor;
%           bound_in(i,2) = bin_y*(1-factor)+bound_in(i,2)*factor;
          

%           figure(1)
%           plot(bound_out(i,1),bound_out(i,2),'b*')
%           plot(bound_in(i,1),bound_in(i,2),'bx')
%           pause
        end        

%  Closes outer loop 
         end; 
%{
figure(2)
clf
hold on
plot(bound_in(:,1),bound_in(:,2),'mx')
plot(bound_out(:,1),bound_out(:,2),'m*')
pause
%}
         for i=1:nbb;
            for j=1:Col;
                  elemnum=(level-1)*nblocks*Col+Col*(block1+i-2)+j;               
                  k1=gl;       
                  for gr=1:3; for gc=1:3;
                     kk=gc+(j-1)*2+Col*2*(i-1);                
                     ex(gr,gc,k1)=bound_in(kk,1)+z(gr)*(bound_out(kk,1)-bound_in(kk,1));
                     ey(gr,gc,k1)=bound_in(kk,2)+z(gr)*(bound_out(kk,2)-bound_in(kk,2));
                  end; end;   
               
                  ex(2,2,k1)=(ex(1,2,k1)+ex(3,2,k1)+ex(2,1,k1)+ex(2,3,k1))*0.5-(ex(1,1,k1)+ex(3,1,k1)+ex(1,3,k1)+ex(3,3,k1))*0.25;
                  ey(2,2,k1)=(ey(1,2,k1)+ey(3,2,k1)+ey(2,1,k1)+ey(2,3,k1))*0.5-(ey(1,1,k1)+ey(3,1,k1)+ey(1,3,k1)+ey(3,3,k1))*0.25;
                
                for gc=1:3;
                  for gr=1:3;
                  pmeshx(elemnum,gr,gc,gl)=ex(gr,gc,gl);
                  pmeshy(elemnum,gr,gc,gl)=ey(gr,gc,gl);
                  pmeshz(elemnum,gr,gc,gl)=zelems(1,gc,gl,1,j,i+block1-1,level);
                  

                  %HERE 2
                  %figure(200)
                  %hold on
                  %plot(squeeze(pmeshx(:,gr,gc,gl)),squeeze(pmeshy(:,gr,gc,gl)),'mx')
                  %pause
                  end;
                end
               
            end
         end
         %
        
         %level
         %gl
         %p
         %pause
         
         

            end;
      %plot(pmeshx(:,:,:,1),pmeshy(:,:,:,1),'mx')
      %pause
end;
end

for level=1:(Laydo-1)
for p=1:npins   
gl=3;

         block1=find(pinmap==p,1,'first');
         blockf=find(pinmap==p,1,'last'); 
         nbb=bb(p);

         for i=1:nbb;
            for j=1:Col;
                  elemnum=(level-1)*nblocks*Col+Col*(block1+i-2)+j;               
                  elemnum1=(level+1-1)*nblocks*Col+Col*(block1+i-2)+j;        
                for gc=1:3;
                  for gr=1:3;
                  pmeshx(elemnum,gr,gc,3)=pmeshx(elemnum1,gr,gc,1);
                  pmeshy(elemnum,gr,gc,3)=pmeshy(elemnum1,gr,gc,1);
                  pmeshz(elemnum,gr,gc,3)=zelems(1,gc,1,1,j,i+block1-1,level+1);    
                  end;
                end
               
            end
         end

end;
end;


level=Laydo;
for p=1:npins   
gl=3;

         block1=find(pinmap==p,1,'first');
         blockf=find(pinmap==p,1,'last'); 
         nbb=bb(p);

         for i=1:nbb;
            for j=1:Col;
                  elemnum=(level-1)*nblocks*Col+Col*(block1+i-2)+j;               
                  elemnum1=Col*(block1+i-2)+j;        
                for gc=1:3;
                  for gr=1:3;
                  pmeshx(elemnum,gr,gc,3)=pmeshx(elemnum1,gr,gc,1);
                  pmeshy(elemnum,gr,gc,3)=pmeshy(elemnum1,gr,gc,1);
                  pmeshz(elemnum,gr,gc,3)=zelems(1,gc,1,1,j,i+block1-1,level)+zelems(1,gc,1,1,j,i+block1-1,2);
                  
                  end;
                end
               
            end
         end

end;


% HERE 3
%{
for ii = 1:elemnum
if (abs(ii- (206938-138240) )<20)
ii
figure(300)
hold on
plot(squeeze(pmeshx(ii,1,:,1)),squeeze(pmeshy(ii,1,:,1)),'mx')
plot(squeeze(pmeshx(ii,2:3,:,1)),squeeze(pmeshy(ii,2:3,:,1)),'kx')
pause
end
end
%}


toc
save savesolid9 pmeshx pmeshy pmeshz r_data

