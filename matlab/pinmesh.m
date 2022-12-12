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
%Df=.2/6.55;
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
%{
track=1;
prev=1;
for i=1:size(pin_blcko,2)
   if pin_blcko(i)>0
      if pin_blcko(i)-prev<=1
         pinmap(track)=pin_blcko(i);
      else
         for j=prev+1:pin_blcko(i)-1
            for k=1:6
               pinmap(track)=j;
               track=track+1;
            end
         end
         pinmap(track)=pin_blcko(i);
      end
      track=track+1;
      prev=pin_blcko(i);
   end
end
pinmap
%}

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

%%%%%%%% HERE LANDON FIX %%%%%%%%%%
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

               %if ~(angle_p<pmax+th && angle_p>pmin+th) && radius<Do/2+1e-2
               if ~(angle_p<pmax+th+pi/24 && angle_p>pmin+th-pi/24) && radius<Do/2+1e-2
                  vcc=1;
               elseif  ~(angle_p<pmax+th && angle_p>pmin+th) && radius<Do/2+1e-2
                   %angle_p-(pmin+th-pi/6)
                   %angle_p-(pmax+th+pi/6)
                  au=min(abs(angle_p-(pmax+th+pi/24)),abs(angle_p-(pmin+th-pi/24)));
                  vcc=1-au/(pi/24)*0.5;
                  
                  if au>pi/12
                      'au is wrong'
                      level
                      p
                     d
                     au
                     %pause
                  end
                      
                  
               elseif radius>Do/2 && radius<Do/2+Dwo/4
                  vcc=0.5-sin(pi/2*(radius-Do/2)/(Dwo/4))*(0.25);
               else
                  vcc=0.25;    
               end;
               if vcc<0.25
                   'vcc is messed up'
                   level
                   p
                   d
                   vcc
                   %pause
               end

                   bdr=sqrt(bdx*bdx+bdy*bdy);
                   bdx=bdx/bdr;
                   bdy=bdy/bdr;

                   %bdx1=bdx*cos(pi*0.5)-bdy*sin(pi*0.5);
                   %bdy1=bdx*sin(pi*0.5)+bdy*cos(pi*0.5);
                   
                   bdx1=-bdy;
                   bdy1=bdx;
                                      
                   bound_in(i,1)=bound_out(i,1)+bdx1*vcc*((Do-Di)/2.0);
                   bound_in(i,2)=bound_out(i,2)+bdy1*vcc*((Do-Di)/2.0);
                   
                   if radius>Do/2+Dw/4 && radius<Do/2+Dw/2
                       theta=atan2(bound_out(i,2)-Ywm,bound_out(i,1)-Xwm);
                       factor=(radius-(Do/2+Dw/4))/(Dw/4);                       
                       bound_ina(i,1)=Xwm+(Dw/2-vcc*(Do-Di)/2.0)*cos(theta);
                       bound_ina(i,2)=Ywm+(Dw/2-vcc*(Do-Di)/2.0)*sin(theta);
                       
                       bound_in(i,1)=factor*bound_ina(i,1)+(1-factor)*bound_in(i,1);
                       bound_in(i,2)=factor*bound_ina(i,2)+(1-factor)*bound_in(i,2);                       
                       
                   end
                   
                   
                   if radius>=Do/2+Dw/2
                       
                       theta=atan2(bound_out(i,2)-Ywm,bound_out(i,1)-Xwm);
                       bound_in(i,1)=Xwm+(Dw/2-vcc*(Do-Di)/2.0)*cos(theta);
                       bound_in(i,2)=Ywm+(Dw/2-vcc*(Do-Di)/2.0)*sin(theta);
                   end
                   
                   
                   if i>1 && i<npts
                       
                       %{
                       hold on
                       plot(bound_in(i,1),bound_in(i,2),'bx')
                       %plot(bound_out(i-1,1),bound_out(i-1,2),'mx')
                       plot(bound_out(i,1),bound_out(i,2),'gx')
                       %plot(bound_out(i+1,1),bound_out(i+1,2),'kx')
                           axis equal
                       pause
                       %}
                     %  
                       
                       % Getting rid of bad angles in fillet
                    rad=((bound_out(i,2)-bound_in(i,2))^2+(bound_out(i,1)-bound_in(i,1))^2)^(1/2);
                    thbef=atan2(bound_out(i-1,2)-bound_out(i,2),bound_out(i-1,1)-bound_out(i,1));
                    thaft=atan2(bound_out(i+1,2)-bound_out(i,2),bound_out(i+1,1)-bound_out(i,1));
                    thi=atan2(bound_in(i,2)-bound_out(i,2),bound_in(i,1)-bound_out(i,1));
                    thbefr=[thbef-2*pi,thbef,thbef+2*pi];
                    thaftr=[thaft-2*pi,thaft,thaft+2*pi];
                    thbef=thbefr(abs(thbefr-thi)<pi);
                    thaft=thaftr(abs(thaftr-thi)<pi);
                    thb=thbef-thi;
                    thc=-(thaft-thi);
                    if thb>3*pi/5
                        %'true1'

                       
                                                     
                        thi=thi+thb-3*pi/5;
                        bound_in(i,1)=bound_out(i,1)+rad*cos(thi);
                        bound_in(i,2)=bound_out(i,2)+rad*sin(thi);
                        %pause
                        

                    end
                    if thc>3*pi/5
                        %'true 2'

                        thi=thi-(thc-3*pi/5);
                        bound_in(i,1)=bound_out(i,1)+rad*cos(thi);
                        bound_in(i,2)=bound_out(i,2)+rad*sin(thi);
                        
                    end

                    %}
                    
                   end

                      
         end
         bound_in2=bound_in;
        
%{
         figure(100) 
%%%%%%%% HERE
         hold on
         plot(bound_in(:,1),bound_in(:,2),'bx')
         plot(bound_out(:,1),bound_out(:,2),'mx')
         plot(bound_in(1,1),bound_in(1,2),'k*')
         plot(bound_out(1,1),bound_out(1,2),'k*')
         pause
         clf
%}

         n_iter=0;
         inege=1.0;
         ibd=zeros(npts);      
       
               
         dxtps=0.001;
         Nopts=(npts-1)*2;
         clear xx0;
         clear df_dx;
         clear ejac;
         
         for i=1:(npts-1);
                 xx0(i)=bound_in(i+1,1);
                 xx0(i+(npts-1))=bound_in(i+1,2);
         end;
         
         xx_i=xx0;
         
         ejac_t=1.0;
         
         %{
         while(ejac_t>0.0) 
                        
         for i=1:(npts-1);
                 bound_in(i+1,1)=xx0(i);
                 bound_in(i+1,2)=xx0(i+(npts-1));
         end;       
         bound_in(1,1)=bound_in(npts,1);    
         bound_in(1,2)=bound_in(npts,2);    
         
         lambda=0.001;    
         n_iter=n_iter+1;
         
         
         if (n_iter==301)
         error('Negative Jacobians, I must stop ...');
         end;
         
         
         bound_in2=bound_in;

         inege=0.0;
         
         for d=1:nelp;             
         for gr=1:3; for gc=1:3; for kk=1:3; 
                     ex(gr,gc,kk)=bound_in(gc+2*(d-1),1)+z(gr)*(bound_out(gc+2*(d-1),1)-bound_in(gc+2*(d-1),1));
                     ey(gr,gc,kk)=bound_in(gc+2*(d-1),2)+z(gr)*(bound_out(gc+2*(d-1),2)-bound_in(gc+2*(d-1),2));
         end; end; end;
 
         for kk=1:3;
         ex(2,2,kk)=(ex(1,2,kk)+ex(3,2,kk)+ex(2,1,kk)+ex(2,3,kk))*0.5-(ex(1,1,kk)+ex(3,1,kk)+ex(1,3,kk)+ex(3,3,kk))*0.25;
         ey(2,2,kk)=(ey(1,2,kk)+ey(3,2,kk)+ey(2,1,kk)+ey(2,3,kk))*0.5-(ey(1,1,kk)+ey(3,1,kk)+ey(1,3,kk)+ey(3,3,kk))*0.25;
         end;
         
         for k=1:3; for j=1:3; for i=1:3; ei(i)=ex(i,j,k); end; dei=Dt*(ei'); for i=1:3; dxdi(i,j,k)=dei(i); end;end;end; 
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ex(i,j,k); end; dej=Dt*(ej'); for j=1:3; dxdj(i,j,k)=dej(j); end;end;end; 
         for k=1:3; for j=1:3; for i=1:3; ei(i)=ey(i,j,k); end; dei=Dt*(ei'); for i=1:3; dydi(i,j,k)=dei(i); end;end;end; 
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ey(i,j,k); end; dej=Dt*(ej'); for j=1:3; dydj(i,j,k)=dej(j); end;end;end;          
         jac=dxdi.*dydj-dxdj.*dydi;                
         for i=1:3; for j=1:3; 
             if (jac(i,j,gl)<1e-7)
                 inege=inege+1.0;  
              end;
    
         end;end;
         ejac(d)=min(min(jac(:,:,gl)));
         end;
         ejac_t=1000.0*sqrt(sum(((ejac-eps_d2).^2).*(ejac<eps_d2)));
         r_jac_p(gl,p,level)=min(ejac);
         ejac_t_o=ejac_t;
         
         
         for k1=1:Nopts             
         xx1=xx0;
         xx1(k1)=xx0(k1)+dxtps;
         for i=1:(npts-1);
                 bound_in(i+1,1)=xx1(i);
                 bound_in(i+1,2)=xx1(i+(npts-1));
         end;    
         bound_in(1,1)=bound_in(npts,1);    
         bound_in(1,2)=bound_in(npts,2);   
         
         for d=1:nelp;             
         for gr=1:3; for gc=1:3; for kk=1:3; 
                     ex(gr,gc,kk)=bound_in(gc+2*(d-1),1)+z(gr)*(bound_out(gc+2*(d-1),1)-bound_in(gc+2*(d-1),1));
                     ey(gr,gc,kk)=bound_in(gc+2*(d-1),2)+z(gr)*(bound_out(gc+2*(d-1),2)-bound_in(gc+2*(d-1),2));
         end; end; end;
 
         for kk=1:3;
         ex(2,2,kk)=(ex(1,2,kk)+ex(3,2,kk)+ex(2,1,kk)+ex(2,3,kk))*0.5-(ex(1,1,kk)+ex(3,1,kk)+ex(1,3,kk)+ex(3,3,kk))*0.25;
         ey(2,2,kk)=(ey(1,2,kk)+ey(3,2,kk)+ey(2,1,kk)+ey(2,3,kk))*0.5-(ey(1,1,kk)+ey(3,1,kk)+ey(1,3,kk)+ey(3,3,kk))*0.25;
         end;
         
         for k=1:3; for j=1:3; for i=1:3; ei(i)=ex(i,j,k); end; dei=Dt*(ei'); for i=1:3; dxdi(i,j,k)=dei(i); end;end;end; 
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ex(i,j,k); end; dej=Dt*(ej'); for j=1:3; dxdj(i,j,k)=dej(j); end;end;end; 
         for k=1:3; for j=1:3; for i=1:3; ei(i)=ey(i,j,k); end; dei=Dt*(ei'); for i=1:3; dydi(i,j,k)=dei(i); end;end;end; 
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ey(i,j,k); end; dej=Dt*(ej'); for j=1:3; dydj(i,j,k)=dej(j); end;end;end;          
         jac=dxdi.*dydj-dxdj.*dydi;                
         ejac(d)=min(min(jac(:,:,gl)));
         
         end;
         ejac_t=1000.0*sqrt(sum(((ejac-eps_d2).^2).*(ejac<eps_d2)));
         df_dx(k1)=(ejac_t-ejac_t_o)/dxtps;
         end;
         
         
         if (ejac_t>0.1);
             
              'negative'
            level
            %gl
            p
             
             %{
                     if (n_iter>1) 
                     dg=df_dx-df_dx_o;
                     ddx=xx0-xx_o;
                     den=(dot(dg,dg));
                     
                     if (abs(den)>0.0) 
                     tlam=dot(dg,ddx)/(dot(dg,dg));
                     else
                     tlam=0.01;
                     end;
                     
                     slam=sign(tlam);
                     lambda=slam*min(abs(tlam),0.1);
                     end;    
                     
                     
                     xx_o=xx0;
                     
%                      if (sum(df_dx.^2)>0.0)
                     xx0=xx_o-lambda*df_dx;
%                      else
%                      for i=1:(npts-2); d=floor((i+1)/2); xx0(i)=xx0(i)+(ejac(d)<0.0)*dxtps;end;   
%                      end;   
                     df_dx_o=df_dx;
                     df_dx=xx0.*0.0; 
             %}
         %end;            
         
%         inege
%  Closes outer loop 
         end; 
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
%
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



%{

for level=1:Lay
for p=1:npins   
gl=2;
         clf 
         hold on
         axis equal

         con1=0;
         con2=0;
         c2w=0;
         w2c=0;

         block1=find(pinmap==p,1,'first');
         blockf=find(pinmap==p,1,'last'); 
         
         nbb=bb(p);
         nelp=nbb*Col;   
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
         
         
         for i=1:nbb;
            for j=1:Col;
                elemnum= (level-1)*nblocks*Col+Col*(block1+i-2)+j;                       
                for gc=1:3;
                  kk=gc+(j-1)*2+Col*2*(i-1);
                  bound_in1(kk,1)=pmeshx(elemnum,1,gc,1);
                  bound_in1(kk,2)=pmeshy(elemnum,1,gc,1);
                  bound_in2(kk,1)=pmeshx(elemnum,1,gc,3);
                  bound_in2(kk,2)=pmeshy(elemnum,1,gc,3);
                  bound_out1(kk,1)=pmeshx(elemnum,3,gc,1);
                  bound_out1(kk,2)=pmeshy(elemnum,3,gc,1);
                  bound_out2(kk,1)=pmeshx(elemnum,3,gc,3);
                  bound_out2(kk,2)=pmeshy(elemnum,3,gc,3); 
                end;
            end;
         end;

         
         
         bo_t=zeros((npts-1)*2);
         for i=1:nbb
            for j=1:Col
                  k=2;
                  kk=k+(j-1)*2+Col*2*(i-1);
                  bo_t(kk)=1.0;
                  bo_t(kk+(npts-1))=1.0;                  
            end
         end         
         
         type=zeros(npts);
       
                 
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
         %plot(Xwm,Ywm,	*')

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

               if ~(angle_p<pmax+th && angle_p>pmin+th) && radius<Do/2+1e-2
               vcc=1.0;
               else
               vcc=0.75;    
               end;
                   
                   bdr=sqrt(bdx*bdx+bdy*bdy);
                   bdx=bdx/bdr;
                   bdy=bdy/bdr;

                   bdx1=bdx*cos(pi*0.5)-bdy*sin(pi*0.5);
                   bdy1=bdx*sin(pi*0.5)+bdy*cos(pi*0.5);

                   bound_in(i,1)=bound_out(i,1)+bdx1*vcc*((Do-Di)/2.0);
                   bound_in(i,2)=bound_out(i,2)+bdy1*vcc*((Do-Di)/2.0);
                      
         end
         bound_in2=bound_in;
                 
         
         n_iter=0;
         inege=1.0;
         ibd=zeros(npts);      
       
         
         dxtps=0.001;
         Nopts=(npts-1)*2;
         clear xx0;
         clear df_dx;
         clear ejac;
         
         for i=1:(npts-1);
                 xx0(i)=bound_in(i+1,1);
                 xx0(i+(npts-1))=bound_in(i+1,2);
         end;
         
         xx_i=xx0;
         
         ejac_t=1.0;
         
         while(ejac_t>0.0) 
             
         for i=1:(npts-1);
                 bound_in(i+1,1)=xx0(i);
                 bound_in(i+1,2)=xx0(i+(npts-1));
         end;       
         
         bound_in(1,1)=bound_in(npts,1);    
         bound_in(1,2)=bound_in(npts,2);
         
         
         for kk=1:(npts-1);
         if (bo_t(kk)>0.9)    
         bound_in(kk,1)=(bound_in1(kk,1)+bound_in2(kk,1)+bound_in(kk+1,1)+bound_in(kk-1,1))*0.5-(bound_in1(kk+1,1)+bound_in1(kk-1,1)+bound_in2(kk+1,1)+bound_in2(kk-1,1))*0.25;
         bound_in(kk,2)=(bound_in1(kk,2)+bound_in2(kk,2)+bound_in(kk+1,2)+bound_in(kk-1,2))*0.5-(bound_in1(kk+1,2)+bound_in1(kk-1,2)+bound_in2(kk+1,2)+bound_in2(kk-1,2))*0.25;
         end;
         end;
                  
         lambda=0.001;    
         n_iter=n_iter+1;
         
         
         if (n_iter==301)
         error('Negative Jacobians, I must stop ...');
         end;
         
         
         bound_in2=bound_in;

         inege=0.0;
         
         for d=1:nelp;             
         for gr=1:3; for gc=1:3; 
                     ex(gr,gc,2)=bound_in(gc+2*(d-1),1)+z(gr)*(bound_out(gc+2*(d-1),1)-bound_in(gc+2*(d-1),1));
                     ey(gr,gc,2)=bound_in(gc+2*(d-1),2)+z(gr)*(bound_out(gc+2*(d-1),2)-bound_in(gc+2*(d-1),2));
         end; end; 
         for gr=1:3; for gc=1:3; 
                     ex(gr,gc,1)=bound_in1(gc+2*(d-1),1)+z(gr)*(bound_out1(gc+2*(d-1),1)-bound_in1(gc+2*(d-1),1));
                     ey(gr,gc,1)=bound_in1(gc+2*(d-1),2)+z(gr)*(bound_out1(gc+2*(d-1),2)-bound_in1(gc+2*(d-1),2));
         end; end; 
         for gr=1:3; for gc=1:3; 
                     ex(gr,gc,3)=bound_in2(gc+2*(d-1),1)+z(gr)*(bound_out2(gc+2*(d-1),1)-bound_in2(gc+2*(d-1),1));
                     ey(gr,gc,3)=bound_in2(gc+2*(d-1),2)+z(gr)*(bound_out2(gc+2*(d-1),2)-bound_in2(gc+2*(d-1),2));
         end; end;      
 
         kk=2;
         ex(2,1,kk)=(ex(1,1,kk)+ex(1,3,kk)+ex(2,1,kk+1)+ex(2,1,kk-1))*0.5-(ex(1,1,kk-1)+ex(3,1,kk-1)+ex(1,1,kk+1)+ex(3,1,kk+1))*0.25;
         ey(2,1,kk)=(ey(1,1,kk)+ey(1,3,kk)+ey(2,1,kk+1)+ey(2,1,kk-1))*0.5-(ey(1,1,kk-1)+ey(3,1,kk-1)+ey(1,1,kk+1)+ey(3,1,kk+1))*0.25; 
         
         kk=2;
         ex(2,3,kk)=(ex(1,3,kk)+ex(3,3,kk)+ex(2,3,kk+1)+ex(2,3,kk-1))*0.5-(ex(1,3,kk-1)+ex(3,3,kk-1)+ex(1,3,kk+1)+ex(3,3,kk+1))*0.25;
         ey(2,3,kk)=(ey(1,3,kk)+ey(3,3,kk)+ey(2,3,kk+1)+ey(2,3,kk-1))*0.5-(ey(1,3,kk-1)+ey(3,3,kk-1)+ey(1,3,kk+1)+ey(3,3,kk+1))*0.25;
 
         kk=2;
         ex(2,2,kk)=(ex(1,2,kk)+ex(3,2,kk)+ex(2,1,kk)+ex(2,3,kk))*0.5-(ex(1,1,kk)+ex(3,1,kk)+ex(1,3,kk)+ex(3,3,kk))*0.25;
         ey(2,2,kk)=(ey(1,2,kk)+ey(3,2,kk)+ey(2,1,kk)+ey(2,3,kk))*0.5-(ey(1,1,kk)+ey(3,1,kk)+ey(1,3,kk)+ey(3,3,kk))*0.25;
         
         
         for k=1:3; for j=1:3; for i=1:3; ei(i)=ex(i,j,k); end; dei=Dt*(ei'); for i=1:3; dxdi(i,j,k)=dei(i); end;end;end; 
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ex(i,j,k); end; dej=Dt*(ej'); for j=1:3; dxdj(i,j,k)=dej(j); end;end;end; 
         for k=1:3; for j=1:3; for i=1:3; ei(i)=ey(i,j,k); end; dei=Dt*(ei'); for i=1:3; dydi(i,j,k)=dei(i); end;end;end; 
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ey(i,j,k); end; dej=Dt*(ej'); for j=1:3; dydj(i,j,k)=dej(j); end;end;end;          
         jac=dxdi.*dydj-dxdj.*dydi;                
         for i=1:3; for j=1:3; 
             if (jac(i,j,gl)<1e-7) 
                 inege=inege+1.0;                
             end;
         end;end;
         ejac(d)=min(min(jac(:,:,gl)));
         end;
         ejac_t=1000.0*sqrt(sum(((ejac-eps_d2).^2).*(ejac<eps_d2)))
         r_jac_p(gl,p,level)=min(ejac);
         ejac_t_o=ejac_t;
         
         
         for k1=1:Nopts             
         xx1=xx0;
         
         if (bo_t(k1)<0.9);
         xx1(k1)=xx0(k1)+dxtps;
         end;
         
         for i=1:(npts-1);
                 bound_in(i+1,1)=xx1(i);
                 bound_in(i+1,2)=xx1(i+(npts-1));
         end;    
         
         bound_in(1,1)=bound_in(npts,1);    
         bound_in(1,2)=bound_in(npts,2); 
         
         for kk=1:(npts-1);
         if (bo_t(kk)>0.9)    
         bound_in(kk,1)=(bound_in1(kk,1)+bound_in2(kk,1)+bound_in(kk+1,1)+bound_in(kk-1,1))*0.5-(bound_in1(kk+1,1)+bound_in1(kk-1,1)+bound_in2(kk+1,1)+bound_in2(kk-1,1))*0.25;
         bound_in(kk,2)=(bound_in1(kk,2)+bound_in2(kk,2)+bound_in(kk+1,2)+bound_in(kk-1,2))*0.5-(bound_in1(kk+1,2)+bound_in1(kk-1,2)+bound_in2(kk+1,2)+bound_in2(kk-1,2))*0.25;
         end;
         end;         
         
         for d=1:nelp;             
         for gr=1:3; for gc=1:3; 
                     ex(gr,gc,2)=bound_in(gc+2*(d-1),1)+z(gr)*(bound_out(gc+2*(d-1),1)-bound_in(gc+2*(d-1),1));
                     ey(gr,gc,2)=bound_in(gc+2*(d-1),2)+z(gr)*(bound_out(gc+2*(d-1),2)-bound_in(gc+2*(d-1),2));
         end; end; 
         for gr=1:3; for gc=1:3; 
                     ex(gr,gc,1)=bound_in1(gc+2*(d-1),1)+z(gr)*(bound_out1(gc+2*(d-1),1)-bound_in1(gc+2*(d-1),1));
                     ey(gr,gc,1)=bound_in1(gc+2*(d-1),2)+z(gr)*(bound_out1(gc+2*(d-1),2)-bound_in1(gc+2*(d-1),2));
         end; end; 
         for gr=1:3; for gc=1:3; 
                     ex(gr,gc,3)=bound_in2(gc+2*(d-1),1)+z(gr)*(bound_out2(gc+2*(d-1),1)-bound_in2(gc+2*(d-1),1));
                     ey(gr,gc,3)=bound_in2(gc+2*(d-1),2)+z(gr)*(bound_out2(gc+2*(d-1),2)-bound_in2(gc+2*(d-1),2));
         end; end;      
 
         kk=2;
         ex(2,1,kk)=(ex(1,1,kk)+ex(1,3,kk)+ex(2,1,kk+1)+ex(2,1,kk-1))*0.5-(ex(1,1,kk-1)+ex(3,1,kk-1)+ex(1,1,kk+1)+ex(3,1,kk+1))*0.25;
         ey(2,1,kk)=(ey(1,1,kk)+ey(1,3,kk)+ey(2,1,kk+1)+ey(2,1,kk-1))*0.5-(ey(1,1,kk-1)+ey(3,1,kk-1)+ey(1,1,kk+1)+ey(3,1,kk+1))*0.25; 
         
         kk=2;
         ex(2,3,kk)=(ex(1,3,kk)+ex(3,3,kk)+ex(2,3,kk+1)+ex(2,3,kk-1))*0.5-(ex(1,3,kk-1)+ex(3,3,kk-1)+ex(1,3,kk+1)+ex(3,3,kk+1))*0.25;
         ey(2,3,kk)=(ey(1,3,kk)+ey(3,3,kk)+ey(2,3,kk+1)+ey(2,3,kk-1))*0.5-(ey(1,3,kk-1)+ey(3,3,kk-1)+ey(1,3,kk+1)+ey(3,3,kk+1))*0.25;
 
         kk=2;
         ex(2,2,kk)=(ex(1,2,kk)+ex(3,2,kk)+ex(2,1,kk)+ex(2,3,kk))*0.5-(ex(1,1,kk)+ex(3,1,kk)+ex(1,3,kk)+ex(3,3,kk))*0.25;
         ey(2,2,kk)=(ey(1,2,kk)+ey(3,2,kk)+ey(2,1,kk)+ey(2,3,kk))*0.5-(ey(1,1,kk)+ey(3,1,kk)+ey(1,3,kk)+ey(3,3,kk))*0.25;
         
         for k=1:3; for j=1:3; for i=1:3; ei(i)=ex(i,j,k); end; dei=Dt*(ei'); for i=1:3; dxdi(i,j,k)=dei(i); end;end;end; 
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ex(i,j,k); end; dej=Dt*(ej'); for j=1:3; dxdj(i,j,k)=dej(j); end;end;end; 
         for k=1:3; for j=1:3; for i=1:3; ei(i)=ey(i,j,k); end; dei=Dt*(ei'); for i=1:3; dydi(i,j,k)=dei(i); end;end;end; 
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ey(i,j,k); end; dej=Dt*(ej'); for j=1:3; dydj(i,j,k)=dej(j); end;end;end;          
         jac=dxdi.*dydj-dxdj.*dydi;                
         ejac(d)=min(min(jac(:,:,gl)));
         end;
         ejac_t=1000.0*sqrt(sum(((ejac-eps_d2).^2).*(ejac<eps_d2)));
         df_dx(k1)=(ejac_t-ejac_t_o)/dxtps;
         end;
         
         
         if (ejac_t>0.0);
             
             'negative' 
             
                     if (n_iter>1) 
                     dg=df_dx-df_dx_o;
                     ddx=xx0-xx_o;
                     den=(dot(dg,dg));
                     
                     if (abs(den)>0.0) 
                     tlam=dot(dg,ddx)/(dot(dg,dg));
                     else
                     tlam=0.01;
                     end;
                     
                     slam=sign(tlam);
                     lambda=slam*min(abs(tlam),0.1);
                     end;    
                     
                     
                     xx_o=xx0;
                     
%                      if (sum(df_dx.^2)>0.0)
                     xx0=xx_o-lambda*df_dx;
%                      else
%                      for i=1:(npts-2); d=floor((i+1)/2); xx0(i)=xx0(i)+(ejac(d)<0.0)*dxtps;end;   
%                      end;   
                     df_dx_o=df_dx;
                     df_dx=xx0.*0.0; 
         end;            
         
%         inege
%  Closes outer loop 
         end; 
       
         for i=1:nbb;
            for j=1:Col;
                  elemnum=(level-1)*nblocks*Col+Col*(block1+i-2)+j;               
                  k1=gl;    
                  
         for gr=1:3; for gc=1:3; 
                     ex(gr,gc,2)=bound_in(gc+2*(j-1),1)+z(gr)*(bound_out(gc+2*(j-1),1)-bound_in(gc+2*(j-1),1));
                     ey(gr,gc,2)=bound_in(gc+2*(j-1),2)+z(gr)*(bound_out(gc+2*(j-1),2)-bound_in(gc+2*(j-1),2));
         end; end; 
         for gr=1:3; for gc=1:3; 
                     ex(gr,gc,1)=bound_in1(gc+2*(j-1),1)+z(gr)*(bound_out1(gc+2*(j-1),1)-bound_in1(gc+2*(j-1),1));
                     ey(gr,gc,1)=bound_in1(gc+2*(j-1),2)+z(gr)*(bound_out1(gc+2*(j-1),2)-bound_in1(gc+2*(j-1),2));
         end; end; 
         for gr=1:3; for gc=1:3; 
                     ex(gr,gc,3)=bound_in2(gc+2*(j-1),1)+z(gr)*(bound_out2(gc+2*(j-1),1)-bound_in2(gc+2*(j-1),1));
                     ey(gr,gc,3)=bound_in2(gc+2*(j-1),2)+z(gr)*(bound_out2(gc+2*(j-1),2)-bound_in2(gc+2*(j-1),2));
         end; end;      
 
         kk=2;
         ex(2,1,kk)=(ex(1,1,kk)+ex(1,3,kk)+ex(2,1,kk+1)+ex(2,1,kk-1))*0.5-(ex(1,1,kk-1)+ex(3,1,kk-1)+ex(1,1,kk+1)+ex(3,1,kk+1))*0.25;
         ey(2,1,kk)=(ey(1,1,kk)+ey(1,3,kk)+ey(2,1,kk+1)+ey(2,1,kk-1))*0.5-(ey(1,1,kk-1)+ey(3,1,kk-1)+ey(1,1,kk+1)+ey(3,1,kk+1))*0.25; 
         
         kk=2;
         ex(2,3,kk)=(ex(1,3,kk)+ex(3,3,kk)+ex(2,3,kk+1)+ex(2,3,kk-1))*0.5-(ex(1,3,kk-1)+ex(3,3,kk-1)+ex(1,3,kk+1)+ex(3,3,kk+1))*0.25;
         ey(2,3,kk)=(ey(1,3,kk)+ey(3,3,kk)+ey(2,3,kk+1)+ey(2,3,kk-1))*0.5-(ey(1,3,kk-1)+ey(3,3,kk-1)+ey(1,3,kk+1)+ey(3,3,kk+1))*0.25;
 
         kk=2;
         ex(2,2,kk)=(ex(1,2,kk)+ex(3,2,kk)+ex(2,1,kk)+ex(2,3,kk))*0.5-(ex(1,1,kk)+ex(3,1,kk)+ex(1,3,kk)+ex(3,3,kk))*0.25;
         ey(2,2,kk)=(ey(1,2,kk)+ey(3,2,kk)+ey(2,1,kk)+ey(2,3,kk))*0.5-(ey(1,1,kk)+ey(3,1,kk)+ey(1,3,kk)+ey(3,3,kk))*0.25;
                
                for gc=1:3;
                  for gr=1:3;
                  pmeshx(elemnum,gr,gc,gl)=ex(gr,gc,gl);
                  pmeshy(elemnum,gr,gc,gl)=ey(gr,gc,gl);
                  pmeshz(elemnum,gr,gc,gl)=zelems(1,gc,gl,1,j,i+block1-1,level);    
                  end;
                end
               
            end
         end
         
        
         level
         p
         

end;
end;
%}


%save savesolid5 pmeshx pmeshy pmeshz


%{


r=1;
size(1:nblocks*Col*R)
size(totelem-nblocks*Col*r+1:totelem)

pmesh_BCs(1:nblocks*Col*R,5)=2;
pmesh_BCs(totelem-nblocks*Col*r+1:totelem,6)=2;
pmesh_BCc(1:nblocks*Col*R)=totelem-nblocks*Col*r+1:totelem;
pmesh_BCc(totelem-nblocks*Col*r+1:totelem)=1:nblocks*Col*R;
toc

size(pmeshx)
size(pmesh_BCs)
size(pmesh_BCi)

for e=1:(nblocks*Col*Laydo)

         for gl=1:3
         for gr=1:3; 
         for gc=1:3; 
                     ex(gr,gc,gl)=pmeshx(e1,gr,gc,gl);
                     ey(gr,gc,gl)=pmeshy(e1,gr,gc,gl);         
 
         end; 
         end; 
         end;
          

         for k=1:3; for j=1:3; for i=1:3; ei(i)=ex(i,j,k); end; dei=Dt*(ei'); for i=1:3; dxdi(i,j,k)=dei(i); end;end;end;
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ex(i,j,k); end; dej=Dt*(ej'); for j=1:3; dxdj(i,j,k)=dej(j); end;end;end;
         for k=1:3; for j=1:3; for i=1:3; ei(i)=ey(i,j,k); end; dei=Dt*(ei'); for i=1:3; dydi(i,j,k)=dei(i); end;end;end;
         for k=1:3; for i=1:3; for j=1:3; ej(j)=ey(i,j,k); end; dej=Dt*(ej'); for j=1:3; dydj(i,j,k)=dej(j); end;end;end;

         jac=dxdi.*dydj-dxdj.*dydi;
         jacvp(elemnum) =min(min(min(jac)));
       
 end; 


save -mat test_psmth2a.mat

fopen( 'pxmesh_simple.out', 'wt' );
fopen( 'pymesh_simple.out', 'wt' );
fopen( 'pzmesh_simple.out', 'wt' );
fopen( 'pxmesh_simple1.out', 'wt' );
fopen( 'pymesh_simple1.out', 'wt' );
fopen( 'pzmesh_simple1.out', 'wt' );


tic
dlmwrite('pxmesh_simple.out',pmeshx,'precision', '%10.6f','-append');
dlmwrite('pymesh_simple.out',pmeshy,'precision', '%10.6f','-append');
dlmwrite('pzmesh_simple.out',pmeshz,'precision', '%10.6f','-append');
dlmwrite('pxmesh_simple1.out',pmeshx,'precision', '%12.8f','-append');
dlmwrite('pymesh_simple1.out',pmeshy,'precision', '%12.8f','-append');
dlmwrite('pzmesh_simple1.out',pmeshz,'precision', '%12.8f','-append');
toc


%}















