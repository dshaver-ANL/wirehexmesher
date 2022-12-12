    function build_wire_template_L(PD,Dw,Df,T)
%
%   Builds an orthogonal grid in the region near a wire of
%   diameter Dw = 0.95 ( P-D ), for D:=1.0.
%
%   Input is PD, the pitch to diameter ratio.
%
%   Typical values:
%
%   PD = 1.135;
%   Dw = 1.00*(PD-1);   %  0.95*(PD-1) is more typical, but here we will 
%   Df = 0.95*(PD-1);   %  "flat" the wire as it nears the neighboring pin.
%



%%                               PD
%%    Current SHARP config:    1.13494;
%%    Todreas & Cheng config:  1.154;
%%    HEDL:                    1.2435    ( = .286/.23)

      close;

      triangulate_wire_template(PD,Dw,Df,T);
%     saves wire_trgl_mesh.mat: t p d wire PD Dw
      %print -deps wire_trgl.eps; display('Printing to wire_trgl.eps')

      %pause

     
      ortho_wire;  
%     loads wire_trgl_mesh.mat: t p d wire PD Dw
%     saves wire_quad_mesh.mat: xx yy tt ss wire_quad_mesh
      %print -deps wire_quad.eps; display('Printing to wire_quad.eps')
      
      %pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function triangulate_wire_template(PD,Dw,Df,TD)

%h = 0.00060 %mesh size
%h = 0.00080 %mesh size
%h = 0.01500; %mesh size %% Use 0.02 for documentation figure
%h = 0.00300; %mesh size
h = 0.00600;

rand('state',111); % Always the same results
set(gcf,'rend','z');

D = 1;                  % Pin diameter
R = D/2;                % Pin radius
P = PD*D;               % Pin pitch
dw = D*Dw;              % small circle for wire
rw = dw/2;              % small circle for wire
rf = Df/2;              % fillet for wire
T = TD*D;

G  = P-2*R;             % Gap between pins
g  = G-dw;              % gap between wire and adjacent pin
Db = 2*R+dw+g/2 + .1*g; % Diameter of big circle
Db = 2*R+5.0*dw;        % Diameter of big circle

% NEW Db
Db=2*(R+dw+1.084/6.55);
Db=2*(R+dw+1.5/6.55);


Rb = Db/2;              % Radius of big circle
Xb = 0.       ;         % Center for big circle



save save1 D R P dw rw rf G g Db Rb

thw= 0;                 % angle of wire
thb= .5*pi;             % domain closing 1/2 angle, opposite wire

fix = zeros(4,2);       % domain fix points

% Plots outline of section
dRb = Xb + (Rb-R); Rt = R  + dRb*(cos(thb)+1)/2.; crct=0;
fix(1,:) = [Rb*cos(thb), Rb*sin(thb)+crct];
fix(2,:) = [ R*cos(thb),  R*sin(thb)];
fix(3,:) = [ R*cos(thb), -R*sin(thb)];
fix(4,:) = [Rb*cos(thb),-Rb*sin(thb)-crct];
plot(fix(:,1),fix(:,2),'ro');axis equal; %pause(.1)




%%%%%%%    GENERATE MESH    %%%%%%%%%%%%
xbox=1.1*(Xb+Db/2);
box = [ -xbox -xbox ; xbox xbox ];  % [ xmin ymin ; xmax ymax ]
wire=[P,R;rw,rf;Db,Xb;thw,thb;fix;T,0];

[p,t]=distmesh2d(@fd1,@fh1,h,box,fix,wire);

post(p,t,@fh1,wire);


d = fd1(p,wire); %% trimesh(t,p(:,1),p(:,2),d); pause

t = chk_jac(t,p);

save wire_trgl_mesh t p d wire PD Dw


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function post(p,t,fh,varargin)

q=simpqual(p,t);
u=uniformity(p,t,fh,varargin{:});
%disp(sprintf(' - Min quality %.2f',min(q)))
%disp(sprintf(' - Uniformity %.1f%%',100*u))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fd1(p,wire) % Signed distance function 

[P,R,rw,rf,Db,Xb,thw,thb,fix]=gwire(wire);

T=wire(9,1);

wire4 = [thw, 2*rw, rf/rw, R, T];

d  =  fdpnwr (p,0.,0.,wire4);
db = -dcircle(p,Xb,0.,Db/2.); % Keep interior of big circle
d  = min(d,db);

%
%  Clip where abs(th) > abs(thb);
%
ct=cos(thb);st=sin(thb);
iyp=find(p(:,2)>0); iym=find(p(:,2)<0);
g=-dline(p,0.,0., ct, st,-9,0); d(iyp)=min(d(iyp),g(iyp));
g=-dline(p,0.,0., ct,-st,-9,0); d(iym)=min(d(iym),g(iym));



%%%%%%%%%%% Landon

%{

pplot=find(d>0);
size(pplot)


skip=1:5:size(pplot);
scatter3(p(pplot(skip),1),p(pplot(skip),2),d(pplot(skip)),5,d(pplot(skip)),'filled')
view(0,90)
axis equal
set(gcf(), 'Renderer', 'painters')
%surf(p(:,1),p(:,2),d)
pause
%}
%%%%%%%%%%%%%%%

d  = -d;   % Per's code likes the negative distance fct.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,R,rw,rf,Db,Xb,thw,thb,fix]=gwire(wire)

P  = wire(1,1);         % Pin pitch
R  = wire(1,2);         % Pin radius

rw = wire(2,1);         % small circle for wire
rf = wire(2,2);         % fillet for wire

Db = wire(3,1);         % Diameter of big circle
Xb = wire(3,2);         % Center for big circle

thw= wire(4,1);         % angle of wire
thb= wire(4,2);         % domain closing 1/2 angle, opposite wire

fix= wire(5:8,:);       % domain fix points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdpnwr(p,xc,yc,wire)
  R = wire(4);  % Pin radius
  d=fdpins(p,xc,yc,R); dw=fdwire(p,xc,yc,wire); d=min(d,dw);
  
  
  %%%%%%%%% Landon
  %dextra=d
  %pasdf = d==dw & d>0;  %Guarantees its outside the pin
 %pasdf = pasdf > -0.001;
  %plot(p(pasdf,1),p(pasdf,2),'r.')
  %pause
  %hold on
  %%%%%%%%%%%%%%%%%%%%%%%%%

  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Standard Wire
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdwire(p,xc,yc,wire)
th  = wire(1);  % Angle of wire
dw  = wire(2);  % wire diameter
rf  = wire(3);  % fillet ratio
rad = wire(4);  % Pin radius
T   = wire(5);  % Wire trim

rw = .5*dw; rf=rf*rw;

xw=xc+(rad+rw)*cos(th); yw=yc+(rad+rw)*sin(th);

%plot(xw,yw,'b*')
%pause
%hold on

thf = ( (rad+rw).^2+(rad+rf).^2-(rw+rf).^2 ) / (2.*(rad+rf)*(rad+rw));
thf = acos(thf);

xtri=zeros(3,1); ytri=zeros(3,1); 




nc=length(xc);
for k=1:nc;
  if k==1; d = dcircle(p,xw(k),yw(k),rw); 
  else;
     g = dcircle(p,xw(k),yw(k),rw); d = min(d,g);
  end;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%% Modified here Landon
  
  line_x=rad+2*rw-T;
  
  these=(p(:,1)>(line_x-T) & abs(p(:,2))<(rw-T)*tan(acos((rw-T)/rw)));
  %these=these(abs(p(:,2))<(rw-T)*tan(acos((rw-T)/rw)));
  d(these)=p(these,1)-line_x;
  %these=(p(:,1)<0.75 & p(:,1)>0.65 & abs(p(:,2))<(rw-T)*tan(acos((rw-T)/rw)));
  %d(these)=p(these


   %plot(0.75,-abs(p(:,2))<(rw-T)*tan(acos((rw-T)/rw)),'b*')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  xtri(1) = xc(k); ytri(1) = yc(k);
  xtri(2) = xw(k); ytri(2) = yw(k);

% Fillets
  xtri(3)=xc(k)+(rad+rf)*cos(th+thf); ytri(3)=yc(k)+(rad+rf)*sin(th+thf);
  dtri=fdedge(p,xtri,ytri); itri=find(dtri > 0);
  n_in_tri=length(itri);
  if n_in_tri > 0;
    df  =-dcircle(p(itri,:),xtri(3),ytri(3),rf); 
    d(itri) = min(d(itri),df);
  end;
  xtri(3)=xc(k)+(rad+rf)*cos(th-thf); ytri(3)=yc(k)+(rad+rf)*sin(th-thf);
  dtri=fdedge(p,xtri,ytri); itri=find(dtri > 0);
  n_in_tri=length(itri);
  if n_in_tri > 0;
    df  =-dcircle(p(itri,:),xtri(3),ytri(3),rf); 
    d(itri) = min(d(itri),df);
  end;

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdpins(p,xc,yc,rad)

    
d = dcircle(p,xc(1),yc(1),rad); nc=length(xc);
for i=2:nc;
   g = dcircle(p,xc(i),yc(i),rad); d = min(d,g);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdedge(p,xe,ye)

    
d=0*p(:,1) + 1.e20; ne=length(xe); 

xc=sum(xe)/ne; yc=sum(ye)/ne;
x1=xe(ne);     y1=ye(ne);
for i=1:ne;
   x0=x1;y0=y1; x1=xe(i); y1=ye(i);
   g = dline(p,x0,y0,x1,y1,xc,yc); d=min(d,g);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=fh1(p,wire) % Size 

   d0  = -fd1(p,wire); d0=max(0,d0);
   e1  = ones(size(p,1),1);

   sh  = 0.10; % shift distance

   for pass=1:4;

     xs  = [sh*e1,.0*e1];
     dx  = -fd1(p+xs,wire); dx=max(0,dx); d0=max(d0,dx);
     dx  = -fd1(p-xs,wire); dx=max(0,dx); d0=max(d0,dx);

     xs  = [.0*e1,sh*e1];
     dx  = -fd1(p+xs,wire); dx=max(0,dx); d0=max(d0,dx);
     dx  = -fd1(p-xs,wire); dx=max(0,dx); d0=max(d0,dx);

     sh = 2*sh/3.;

   end;

   alpha = 0.020;
   h = min(d0,alpha)/alpha;

   th = atan2(p(:,2),p(:,1))*(2./pi); th = abs(th);
   th = (th/.3).^2;  th = min(th,1); th = (th+.5)/(1.5);
   h = min(h,th);

   h = .18 + h;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ortho_wire;

set(gcf,'rend','z');

hold off;
load wire_trgl_mesh; %% t p d wire PD Dw

po = p; to=t; nb=size(po,1);

[P,R,rw,rf,Db,Xb,thw,thb,fix]=gwire(wire);

ng = size(fix,1);  % # boundary segments, Gamma


dg = zeros(ng,2);  % flag for Dirichlet boundaries and bdry values
dg(1,1) = 1; dg(3,1) = 1;  % Gamma 1 & 3 are Dirichlet
dg(1,2) = 1; dg(3,2) = 0;  % u(g1) = 0;  u(g3) = 1;
source=0; % No source for first laplace solve
[ub,edge,efix]=laplace_fem(to,po,fix,dg,source,R,rw);


gamma4=edge(efix(4):efix(5));

th=atan2(po(gamma4,2),po(gamma4,1));u4=ub(gamma4); 
%hold off;plot(th,u4)

%
%  Spline
%

%m=991; % RESOLUTION
m=9991;


t0=min(th); t1=max(th); dt=(t1-t0)/m; xx=dt*(0:m)+t0; 
xm=.5*(xx(1:m)+xx(2:m+1));Dt=t1-t0; dt=dt*(Dt+3*abs(xm).^5+0.025./(.08+xm.^4)); 
xx=[0,cumsum(dt)]; X0=min(xx);X1=max(xx);
xx = t0+(t1-t0)*(xx-X0)/(X1-X0);


uu=spline(th,u4,xx); uu(1)=1.e-7;uu(m+1)=1-1.e-7; tt=xx;


dg = zeros(ng,2);  % flag for Dirichlet boundaries and bdry values
dg(2,1) = 1; dg(4,1) = 1;  % Gamma 2 & 4 are Dirichlet
dg(2,2) = 1; dg(4,2) = 0;  % u(g2) = 1;  u(g4) = 0;
source=0; % Include source or not
[u2,edge,efix]=laplace_fem(to,po,fix,dg,source,R,rw);

%%
%% Uncomment below to see a mesh;
%%

% hold on; 
% for k=1:m+1; a=uu(k); x1=cont_tri(to,po,ub,a); axis equal; end;
% for a=.2:.2:.8; x2=cont_tri(to,po,u2,a); axis equal; end;
% pause % pause(2)

%ms = 8;
ms = 8; xo=[po,u2]; ssl=[0:ms]'/ms; xx=zeros(ms+1,m+1);yy=xx;
for k=1:m+1; 
   a=uu(k); x1=cont_tri_val(to,xo,ub,a); 
   x1=vector_clean(x1,3);
   x = x1(:,1); y = x1(:,2); s = x1(:,3);
   
   x = smooth(x1(:,1),2,.3); y = smooth(x1(:,2),2,.3); s =smooth(x1(:,3),2,.3); 

   %s=(0:1/(size(s,1)-1):1)';

   xx(:,k) = spline(s,x,ssl);  xx(:,k) = smooth(xx(:,k),1,.5);
   yy(:,k) = spline(s,y,ssl);  yy(:,k) = smooth(yy(:,k),1,.5);
end;

% Get s distribution that results in evenly spaced circumferential lines
ymin=min(yy(:,1));
ymax=max(yy(:,1));
ydist=(ymin+(ymax-ymin)/(ms)*(0:ms))';
snew=spline(yy(:,1),ssl,ydist);
snew2=flipud(-snew+1); % This flips it
%ssl=snew;
ssl=snew2;


x=xx; y=yy; t=tt; s=ssl; close all; hold on; axis equal; axis off;


% Distribute points across inside boundary 
% and find theta of corresponding outside position
% this is the theta of the inside point

ratio=1.03; % Ratio for geometric series shifting inside distribution
power=.75; % Power used to weight the combination of inside and outside when merging meshes
Lsmoothe=1; % Number of Laplacian Smoothing routines

nt=100; t0=min(t);t1=max(t); % nt = NUMBER OF LINES
ns=length(s); xt=zeros(ns,nt+1);yt=xt;

%   Shift distribution toward fillets

% Determine location of fillet
A = R+rw;
B = R+rf;
C = rw+rf;

aa = acos((B^2 + C^2 - A^2)/(2*B*C));
ab = acos((A^2 + C^2 - B^2)/(2*A*C));
ac = acos((A^2 + B^2 - C^2)/(2*A*B));

xth=B*cos(-ac)+rf/2*cos(ab+aa/2);
yth=B*sin(-ac)+rf/2*sin(ab+aa/2);


% Find distance from edge to fillet and fillet to center

distf=0;
distc=0;
dtot=0;
for i=1:m
    d=((x(ns,i+1)-x(ns,i))^2+(y(ns,i+1)-y(ns,i))^2)^(1/2);
    dtot=dtot+d;
    if x(ns,i)<xth && y(ns,i)<0;
        distf=distf+d;
    end
    if x(ns,i)>xth && y(ns,i)<0
        distc=distc+d;
    end
end

% Distribute points concentrated more toward the fillet
dt=zeros(nt,1);
nf=floor(nt*distf/dtot);
nc=ceil(nt*distc/dtot);

r=ratio;
Af=distf*(1-r)/(1-r^(nf));
for i=1:nf
    dt(i)=Af*r^(nf-i);
    dt(nt+1-i)=dt(i);
end
r=(r-1)*2+1; % This should make the wire spacing nearly match up with the pin
Ac=distc*(1-r)/(1-r^(nc));
for i=nf+1:nf+nc
    dt(i)=Ac*r^(i-nf-1);
    dt(nt+1-i)=dt(i);
end

% For each point on the inside, find corresponding theta on outside
xf = zeros(nt,1);
yf=xf;
tt = zeros(1,nt+1);
dleft=0;
j=0;
tt(1)=t0;
tt(nt+1)=t1;
for i=2:nt;
    d=dt(i-1);
    if d>dleft;
        j=j+1;
        d=d-dleft;
        dnext = dxy(x(ns,j),y(ns,j),x(ns,j+1),y(ns,j+1));
        while d>dnext
            j=j+1;
            d=d-dnext;   
            dnext = dxy(x(ns,j),y(ns,j),x(ns,j+1),y(ns,j+1));
        end
        tt(i) = atan(y(1,j)/abs(x(1,j)))+d/dnext*(atan(y(1,j+1)/abs(x(1,j+1)))-atan(y(1,j)/abs(x(1,j))));
        dleft = dnext-d;
    else
        tt(i) = tt(i-1)+d/dleft*(atan(y(1,j+1)/abs(x(1,j+1)))-atan(y(1,j)/abs(x(1,j))));
        dleft=dleft-d;
    end
end

% Uniform outside distribution

t0=min(t);t1=max(t); dt=(t1-t0)/nt; tt2 = t0 + dt*(0:nt); % nt = NUMBER OF LINES
ns=length(s); xt=zeros(ns,nt+1);yt=xt;

% Combine inside and outside distributions to make nicer mesh
xt=zeros(ns,nt+1);yt=xt;xt1=xt;yt1=xt;xt2=xt;yt2=xt;

p=power;

for i=1:ns
    xt1(i,:)=spline(t,x(i,:),tt);
    yt1(i,:)=spline(t,y(i,:),tt); 
    xt2(i,:)=spline(t,x(i,:),tt2);
    yt2(i,:)=spline(t,y(i,:),tt2);
    %xt(i,:)=xt1(i,:)*(i-1)/(ns-1)+xt2(i,:)*(1-(i-1)/(ns-1)); No weight
    %yt(i,:)=yt1(i,:)*(i-1)/(ns-1)+yt2(i,:)*(1-(i-1)/(ns-1));
    xt(i,:)=xt1(i,:)*((i-1)/(ns-1))^(p)+xt2(i,:)*(1-((i-1)/(ns-1))^(p));
    yt(i,:)=yt1(i,:)*((i-1)/(ns-1))^(p)+yt2(i,:)*(1-((i-1)/(ns-1))^(p));
end

%Laplacian smoothing
%{
L=Lsmoothe;
for k=1:L
for i=2:ns-1
    for j=2:nt
        xt(i,j)=( xt(i-1,j)+xt(i+1,j)+xt(i,j-1)+xt(i,j+1) )/4;
        yt(i,j)=( yt(i-1,j)+yt(i+1,j)+yt(i,j-1)+yt(i,j+1) )/4;
    end
end
end
%}

xx=xt;
yy=yt;


wire_quad_mesh = wire;
save wire_quad_mesh_flat   xx yy tt ssl wire_quad_mesh PD Dw
display('Data saved to wire_quad_mesh_flat.mat')

% Display mesh


%{

%for i=1:ns; xt(i,:)=spline(t,x(i,:),tt2); yt(i,:)=spline(t,y(i,:),tt2);end;

ns=50; s0=min(s);s1=max(s); ds=(s1-s0)/ns; ssl = s0 + ds*(0:ns); % ns = LINE RESOLUTION
for j=1:nt+1;xs=spline(s,xt(:,j),ssl);ys=spline(s,yt(:,j),ssl);plot(ys,xs,'r-');end;

%
%   ns lines at constant r;
%
%ns=15;
ns=15; s0=min(s);s1=max(s); ds=(s1-s0)/ns; ssl = s0 + ds*(0:ns); % ns = NUMBER OF LINES
nt=length(t); xs=zeros(ns+1,nt); ys=xs;
for j=1:nt; xs(:,j)=spline(s,x(:,j),ssl); ys(:,j)=spline(s,y(:,j),ssl); end;

nt=999; t0=min(t);t1=max(t); dt=(t1-t0)/nt; tt = t0 + dt*(0:nt);% nt = LINE RESOLUTION
for i=1:ns+1;xt=spline(t,xs(i,:),tt);yt=spline(t,ys(i,:),tt);plot(yt,xt,'r-');end;

pause(4)

%}

hold off
mesh(xx,yy,0*xx,'EdgeColor','blue'); view(0,90); axis image; axis off;
set(gcf(), 'Renderer', 'painters')
pause




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  END OF TEST CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iskin = skin(ti); 
%
%  ordered skin for a list of triangles
%

tmin=min(min(ti));
t  = ti-tmin+1;               % triangle vertices in [1,...,tmax];

tmax=max(max(t));
nt = size(t,1);               % # triangles

a = spalloc(tmax,tmax,3*nt);  % Edge graph


for k=1:3; k1=k+1; if k1>3; k1=1; end;
    v0=min(t(:,k),t(:,k1));
    v1=max(t(:,k),t(:,k1));
    a = a + sparse(v0,v1,1,tmax,tmax);
end;

%spy(a)

[r,c] = find(a==1);           % Faces have only one entry in graph

nr=max(r); nc=max(c); n=max(nr,nc);
nlinks=zeros(n,1); links=zeros(n,2); nedge=length(r); 

for k=1:nedge;
    i=c(k); j=r(k); nlinks(i)=nlinks(i)+1; links(i,nlinks(i))=j;
    i=r(k); j=c(k); nlinks(i)=nlinks(i)+1; links(i,nlinks(i))=j;
end;

iskin=zeros(nedge+1,1);
i=r(1); iskin(1)=i; inext = links(i,1);
for k=2:nedge;
    ilast=i; i=inext; iskin(k)=i;
    inext = links(i,1); if inext==ilast; inext = links(i,2); end;
end;
iskin(nedge+1)=inext; iskin=iskin+tmin-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ub,edge,efix]=laplace_fem(to,po,fix,dg,source,rp,rw);
%
%  Evaluate solution to Laplace's equation
%
%

sstrength=0.08;

p=po;t=to; nb=size(po,1); ng=size(fix,1);  % # boundary segments, Gamma

edge = skin(to); nedge=length(edge)-1;
plot(p(edge,1),p(edge,2),'r-'); axis equal; axis off; hold on; 
xmin=min(p(:,1)); xmax=max(p(:,1));
ymin=min(p(:,2)); ymax=max(p(:,2));
xmin=min(xmin,ymin); xmax=max(xmax,ymax);
axis ([xmin xmax xmin xmax]);

%
%   Enforce fix points to be members of edge set
%
nfix=size(fix,1); e1=[1;1]; ifix = zeros(nfix,1);
for k=1:nfix;
   pf=p(edge,:); 
   pf(:,1)=pf(:,1)-fix(k,1); pf(:,2)=pf(:,2)-fix(k,2); pf=abs(pf); pf=pf*e1;
   [fmin,jfix]=min(pf); ifix(k)=edge(jfix); efix(k)=jfix;
   fix(k,1)=p(ifix(k),1); fix(k,2)=p(ifix(k),2);
end;

%
%  Perform cyclic shift to get ifix(1) in first edge position
%
i1=find(edge==ifix(1));
ec=zeros(nedge+1,1);ec(1:nedge-i1+1)=edge(i1:nedge);
ec(nedge-i1+2:nedge)=edge(1:i1-1);edge=ec;edge(nedge+1)=edge(1);
%
%  If area < 0, flip
%
area = 0;
for k=2:nedge+1;
  i0=edge(k-1);i1=edge(k); dx=p(i1,1)-p(i0,1); area = area+p(i1,2)*dx; %y*dx
end; area=-area; if area < 0; edge(nedge+1:-1:1)=edge; end;

for k=1:nfix;   % Re-identify fix points
   pf=p(edge,:); 
   pf(:,1)=pf(:,1)-fix(k,1); pf(:,2)=pf(:,2)-fix(k,2); pf=abs(pf); pf=pf*e1;
   [fmin,jfix]=min(pf); ifix(k)=edge(jfix); efix(k)=jfix;
   fix(k,1)=p(ifix(k),1); fix(k,2)=p(ifix(k),2);
end; efix(nfix+1)=nedge+1;



%
%  Set Dirichlet BCs on specified Gammas
%


ndbc=0; ub = zeros(nb,1);  % Set up RHS and identify Dirichlet bdrys
for k=1:ng; if dg(k,1)==1; ndbc=ndbc+1;
    if ndbc==1; edbc=edge(efix(k):efix(k+1)); else;
       edbc=union(edbc,edge(efix(k):efix(k+1)));  end;
    ub(edge(efix(k):efix(k+1))) = dg(k,2);
end; end;
nbc=length(edbc); n=nb-nbc;

[t,p,ip] = reorder_p(to,po,edbc); ub=ub(ip); % Dirichlet points ordered last


[Ab,Q]=afem(p,t); % Build 2D FEM stiffness matrix

A  = Ab(1:n,1:n);
R  = speye(nb); R = R(1:n,:);

b  = -R*Ab*ub;

% Force circumferential lines inward

if source==1
    for k=1:n

        if  (p(k,2)<(rw-rp)/(rw+rp)*p(k,1)+rp && p(k,2)>0.1) || ...
            (p(k,2)>(rp-rw)/(rw+rp)*p(k,1)-rp && p(k,2)<-0.1)
            b(k)=b(k)-sstrength;
            plot(p(k,1),p(k,2))
        end
     
    end  
end


u  = A\b;
ub = ub + R'*u;

%trisurf(t,p(:,1),p(:,2),ub); axis equal; drawnow; 


ub(ip) = ub; % Transform back to original ordering

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Contour, plus evaluate another variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,W] = cont_tri2(t,p,w,u,uv); 

minu = min(u(t)');maxu=max(u(t)');

ic = find(minu <= uv & uv < maxu);

hold on;
nt = length(ic); X=zeros(2*nt,2);W=zeros(2*nt); k=1;
for e=1:nt; k1=k+1;
  [x,wt,a] = cont_e2(t,p,w,ic(e),u,uv);
  plot(x(:,1),x(:,2),'r-'); %pause(.1)
  W(k:k1)=wt; X(k:k1,:)=x; k=k+2;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Single triangle contour, plus evaluate another variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,wt,a] = cont_e2(t,X,w,e,u,uv); % single element contour line

  s = u(t(e,:)) ;
  p = X(t(e,:),:);

  k=0;
  a=0;
  for j=1:3;
     j1 = j+1;  if j1>3 , j1=1; end;
     if (s(j)-uv)*(s(j1)-uv) <=0 & s(j) ~= s(j1)
        k = k+1;
        a(k)   = (uv-s(j))/(s(j1)-s(j));
        x(k,:) = p(j,:) + a(k)*(p(j1,:)-p(j,:));
        wt(k)  = w(j)   + a(k)*(w(j1)  -w(j));
     end;
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XX] = cont_tri_val(t,x,u,uv); 
minu = min(u(t)');maxu=max(u(t)'); ic = find(minu <= uv & maxu > uv);
nt=length(ic); ncol=size(x,2); XX=zeros(2*nt,ncol); k=1;
for e=1:nt; [XX(k:k+1,:),a] = cont_e(t,x,ic(e),u,uv); k=k+2; end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xc] = vector_clean(x,c);  % sort by col c, eliminate near duplicates.
if size(x,1) > 0;
  [s,I]=sort(x(:,c)); xc=x(I,:); n=length(I); epsilon = 100*eps; k=1; 
  for i=2:n;if abs(xc(i,c)-xc(k,c))>epsilon;k=k+1;xc(k,:)=xc(i,:);end;end;
  xc=xc(1:k,:);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = smooth(x,npass,a);  % Laplacian smoothing
s=x; n1=size(s,1); n2=size(s,2); if n2 > n1; s=x'; end; n=size(s,1);
s1=s(1); sn=s(n);
e = ones(n,1); A = spdiags([.5*a*e (1-a)*e .5*a*e], -1:1, n, n);
for k=1:npass; s=A*s; s(1)=s1; s(n)=sn; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = cont_tri(t,p,u,uv); 

minu = min(u(t)');maxu=max(u(t)');

ic = find(minu <= uv & uv < maxu);
% hold on;
% plot(p(t(ic,:),1),p(t(ic,:),2),'ro')
% pause;
% x=0;

  nt = length(ic);
  for e=1:nt
    [x,a] = cont_e(t,p,ic(e),u,uv);
    hold on;
%   plot(x(:,1),x(:,2),'g.'); %pause(.1)
    plot(x(:,1),x(:,2),'r-'); %pause(.1)
  end;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [d] = dxy(x1,y1,x2,y2);
     
 d=((x2-x1)^2+(y2-y1)^2)^(1/2);
        

     
        
        
        
        
        
        
        
        
        
        
        
        
        
        