%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xyq=uniform_blocks(xyblock,th,X,Y,pin_block,pin_blcko,type_block);

fprintf('used uniform blocks')
pause

PD = .286/.23;
dw = 0.95*(PD-1);
rf = 1.0;
rad= 0.5;

wire(1)=th;  % Angle of wire
wire(2)=dw;  % wire diameter
wire(3)=rf;  % fillet ratio
wire(4)=rad; % Pin radius


xyq=xyblock;


plot(X,Y,'x')
pause



b0=1;
for pin=1:13;
  b1=b0; while pin_block(b1)==pin_block(b0); b1=b1+1; end; b1=b1-1;
  pin = pin_block(b0);
  pno = pin_blcko(b0);
  xyq(:,:,b0:b1,:)=uni_block(xyblock,type_block,b0,b1,X(pno),Y(pno),wire);
  pin
  b0=b1+1;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xyq=uni_block(xyb,type,b0,b1,XP,YP,wire)
%
%  Build new block file that is nearly uniformly spaced around
%  pin circumference.

%  To minimize distortion, fixed points in the equi-arclength
%  repartition are taken to be the tip of the wire and the 
%  point that is 180 degrees opposite this point.
%
%
%  Equi-length partitions are based on assumed values of ne_r,
%  here taken as: ne_r(1)=8; ne_r(2)=13;ne_r(3)=10;
%  which was found to give a reasonable distribution. 
%
%  Note that the output blocks are of fixed dimension [m,n]
%  and not fixed block-to-block density, although intra-block
%  density is nearly uniform.   build_ess.m relies on uniform
%  intrablock density when partitioning each block into nels x ne_r
%  spectral elements.
%

xyq=xyb(:,:,b0:b1,:);

nb = 1 + (b1-b0);
m=size(xyb,1);
n=size(xyb,2);


ntot = (b1-b0+1)*m*n;
xy0(1) = sum(sum(sum(xyb(:,:,b0:b1,1))))/ntot;  % center of
xy0(2) = sum(sum(sum(xyb(:,:,b0:b1,2))))/ntot;  % gravity


hold off;

X=zeros(m,6*8*12);
Y=zeros(m,6*8*12);
l = 0; ll= 1;
for b=b0:b1;
  l = l+1;
  xy1(1)=xyb(1,1,b,1);
  xy1(2)=xyb(1,1,b,2);
  xym(1)=xyb(1,n,b,1);
  xym(2)=xyb(1,n,b,2);
  [xx,yy,nn,ne_r] = qspl(xyb(:,:,b,1),xyb(:,:,b,2),type(b));
  nn=nn-1;
  nnl(l) = nn;
  X(:,ll:ll+nn-1) = xx(:,1:nn);
  Y(:,ll:ll+nn-1) = yy(:,1:nn);
  if orient(xy0,xy1,xym)==1 ;%% [ l 'flip']
     X(:,ll:ll+nn-1) = xx(:,nn+1:-1:2);
     Y(:,ll:ll+nn-1) = yy(:,nn+1:-1:2);
  end;
  ll = ll+nn;
end;

ll=ll-1; X=X(:,1:ll); Y=Y(:,1:ll);


m = size(X,1);   plot(X(m,:),Y(m,:),'ro-'); axis equal; 
pause %

p=[X(m,:)' , Y(m,:)']; p=arc_shift(p,XP,YP,wire); p=[ p ; p(1,:) ];

l = 0; l0=0; xx=zeros(2,1);yy=xx;

m = size(xyb,1); n = size(xyb,2); 
for b=b0:b1;
  l = l+1;
  xy1(1)=xyb(1,1,b,1);
  xy1(2)=xyb(1,1,b,2);
  xym(1)=xyb(1,n,b,1);
  xym(2)=xyb(1,n,b,2);

  nn=nnl(l);
  l1=l0+1;
  l2=l0+nn+1;
  [xb,yb,sb]=spline_block_arc2(p(l1:l2,1),p(l1:l2,2),n);
  

%
% Must pick up overlap!
%

  if orient(xy0,xy1,xym)==1 ; xb=xb(n:-1:1); yb=yb(n:-1:1); end;

  xyq(m,:,l,1)=xb; xyq(m,:,l,2)=yb;
return
  for j=1:n;
     xx(1)=xyq(1,j,l,1); xx(2)=xyq(m,j,l,1);
     yy(1)=xyq(1,j,l,2); yy(2)=xyq(m,j,l,2);
     [xyproj,grd] = pfdwire([xx(2),yy(2)],XP,YP,wire); % pin projection
     xx(2)=xyproj(1); yy(2)=xyproj(2);
     [xb,yb,sb]=spline_block_arc2(xx,yy,m);
     xyq(:,j,l,1)=xb;
     xyq(:,j,l,2)=yb;%print
  end;

  l0=l2-1;

% if b<11;
%   spline_block(xyq(:,:,l,1),xyq(:,:,l,2)); axis image;%pause(0.01)
% end;

end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
function spline_block(X,Y);

    m=size(X,1); 
    n=size(X,2);
    s=0.1;



%   xx=spline(1:n,X,1:s:n); yy=spline(1:n,Y,1:s:n);
    xx=pchip (1:n,X,1:s:n); yy=pchip (1:n,Y,1:s:n);
    for i=1:m; plot(xx(i,:),yy(i,:),'m-'); hold on; end;
    for i=1:1; plot(xx(i,:),yy(i,:),'k-'); end;

%   xx=spline(1:m,X',1:s:m); yy=spline(1:m,Y',1:s:m);
    xx=pchip (1:m,X',1:s:m); yy=pchip (1:m,Y',1:s:m);

%   for j=1:8:n; plot(xx(j,:),yy(j,:),'m-'); end; pause (.01);
    for j=1:n; plot(xx(j,:),yy(j,:),'m-'); end; pause (.01);

%   for j=1:1; plot(xx(j,:),yy(j,:),'k-'); end; pause (.01);
%   for j=1:1; plot(xx(j,:),yy(j,:),'k.'); end; pause (.01);
    for j=4:4; plot(xx(j,:),yy(j,:),'r.'); end; pause (.01);
    for j=4:4; plot(xx(j,:),yy(j,:),'r-'); end; pause (.01);

    axis image

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = orient(x0,x1,x2);

 o = 0; % right handed

 v1=x1-x0; v2=x2-x0;

 cross = v1(1)*v2(2) - v1(2)*v2(1);

  if cross < 0; o=1; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,yy,nnr,ne_r]=qspl(X,Y,type);

n=size(X,2); s=(0:(n-1))/(n-1);

%ne_r=zeros(3,1); ne_r(1)=8; ne_r(2)=11;ne_r(3)=10; % # elements in azimuth
 ne_r=zeros(3,1); ne_r(1)=8; ne_r(2)=13;ne_r(3)=10; % # elements in azimuth

nnr=8*ne_r(type)+1; ss=(0:(nnr-1))/(nnr-1);



plot(X,Y,'kx'); axis equal; hold on
% pause % see individual sections
 xx=spline(s,X,ss);
 yy=spline(s,Y,ss);
 
%xx=pchip (s,X,ss); yy=pchip (s,Y,ss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xp,dx]=pfdwire(p,xc,yc,wire)
%
%  xp and yp are the projections of p onto the wire.
%
%  dx and dy are the normalized gradients of fdwire, at (xp,yp)
%
%
%

xp=p; dx=p;
n=size(p,1);
for k=1:n
  xp(k,:)=pin_proj(p(k,:),xc,yc,wire);
  dx(k,:)=pin_grad(xp(k,:),xc,yc,wire);
  dnorm=dx(k,1)*dx(k,1)+dx(k,2)*dx(k,2); dx(k,:)=dx(k,:)/sqrt(dnorm);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdwire(p,xc,yc,wire)
th  = wire(1);  % Angle of wire
dw  = wire(2);  % wire diameter
rf  = wire(3);  % fillet ratio
rad = wire(4);  % Pin radius

rw = .5*dw; rf=rf*rw;

xw=xc+(rad+rw)*cos(th); yw=yc+(rad+rw)*sin(th);

thf = ( (rad+rw).^2+(rad+rf).^2-(rw+rf).^2 ) / (2.*(rad+rf)*(rad+rw));
thf = acos(thf);

xtri=zeros(3,1); ytri=zeros(3,1); 

nc=length(xc);
for k=1:nc;
  if k==1; d = dcircle(p,xw(k),yw(k),rw); else;
     g = dcircle(p,xw(k),yw(k),rw); d = min(d,g);
  end;
  
  xtri(1) = xc(k); ytri(1) = yc(k);
  xtri(2) = xw(k); ytri(2) = yw(k);

% Fillets
  xtri(3)=xc(k)+(rad+rf)*cos(th+thf);
ytri(3)=yc(k)+(rad+rf)*sin(th+thf);
  dtri=fdedge(p,xtri,ytri); itri=find(dtri > 0);
  n_in_tri=length(itri);
  if n_in_tri > 0;
    df  =-dcircle(p(itri,:),xtri(3),ytri(3),rf); 
    d(itri) = min(d(itri),df);
  end;
  xtri(3)=xc(k)+(rad+rf)*cos(th-thf);
ytri(3)=yc(k)+(rad+rf)*sin(th-thf);
  dtri=fdedge(p,xtri,ytri); itri=find(dtri > 0);
  n_in_tri=length(itri);
  if n_in_tri > 0;
    df  =-dcircle(p(itri,:),xtri(3),ytri(3),rf); 
    d(itri) = min(d(itri),df);
  end;

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdedge(p,xe,ye)
d=0*p(:,1) + 1.e20; ne=length(xe); 

xc=sum(xe)/ne; yc=sum(ye)/ne;
x1=xe(ne);     y1=ye(ne);
for i=1:ne;
   x0=x1;y0=y1; x1=xe(i); y1=ye(i);
   g = dline(p,x0,y0,x1,y1,xc,yc); d=min(d,g);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=dline(p,x1,y1,x2,y2,xc,yc)

vx = x2-x1; cx = xc-x1;
vy = y2-y1; cy = yc-y1;

vxc = vx*cy-vy*cx;

nx = -vxc*vy;
ny =  vxc*vx;

nn = sqrt(nx*nx+ny*ny);
nx = nx/nn;
ny = ny/nn;

d = nx*(p(:,1)-x1) + ny*(p(:,2)-y1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=dcircle(p,xc,yc,r)

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.

d=sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2)-r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grd=pin_grad(x,xc,yc,wire);
h=1.e-5;
px = [ x(1)-h , x(1)+h ,  x(1)   , x(1)   ]'; 
py = [ x(2)   , x(2)   ,  x(2)-h , x(2)+h ]'; 
X = [ px , py ]; Z = fdpnwr(X,xc,yc,wire);
grd=0*x; grd(1)=.5*(Z(2)-Z(1))/h; grd(2)=.5*(Z(4)-Z(3))/h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdpnwr(p,xc,yc,wire)
  R = wire(4);  % Pin radius
  d=fdpins(p,xc,yc,R); dw=fdwire(p,xc,yc,wire); d=min(d,dw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdpins(p,xc,yc,rad)
d = dcircle(p,xc(1),yc(1),rad); nc=length(xc);
for i=2:nc;
   g = dcircle(p,xc(i),yc(i),rad); d = min(d,g);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=pin_proj(x0,X,Y,wire);

%     Project x0 onto pin+wire

th    = wire(1);  % Angle of wire
dw    = wire(2);  % wire diameter
rf    = wire(3);  % fillet ratio
R     = wire(4);  % Pin radius


x1=x0;
h=fdpnwr(x1,X,Y,wire); grd=pin_grad(x1,X,Y,wire);
iter=0;
while abs(h) > 100*eps & iter < 30; iter=iter+1;
   x0=x1; x1=x0-h*grd/norm(grd); h=fdpnwr(x1,X,Y,wire);
%  [x1, iter, h]
end;
x=x1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q=arc_shift(p,XP,YP,wire);

th    = wire(1);  % Angle of wire
dw    = wire(2);  % wire diameter
rf    = wire(3);  % fillet ratio
R     = wire(4);  % Pin radius

% (XP,YP) = pin coordinates

% p = current set of points on pin boundary
% q = shifted set of points on pin boundary

m=size(p,1);
p0(:,1)=p(:,1)-XP; p0(:,2)=p(:,2)-YP;

th0 = atan2(p0(:,2),p0(:,1)) - th + 6*pi;
for k=1:4; if max(th0)*min(th0) > 0; th0=th0-2*pi; end; end;

[thm,i0]=min(abs(th0));
thm = th0(i0);
if thm<0; i0=i0+1; end; if i0>m; i0=1; end; thm=th0(i0);
ps=cyc_shift(p0,i0,1);

th1 = atan2(-ps(:,2),-ps(:,1)) - th + 6*pi;
for k=1:4; if max(th1)*min(th1) > 0; th1=th1-2*pi; end; end;
[thm,i1]=min(abs(th1)); if thm>0; i1=i1-1; end;

xa=zeros(m,1);ya=zeros(m,1);sa=zeros(m,1);

[xa(1:i1),ya(1:i1),sa(1:i1)]=spline_block_arc1(ps(1:i1,1),ps(1:i1,2)); % ms = 8
i2=i1+1;
[xa(i2:m),ya(i2:m),sa(i2:m)]=spline_block_arc1(ps(i2:m,1),ps(i2:m,2));

q=[ xa+XP , ya+YP ];
q=cyc_shift(q,i0,-1);
%plot(q(:,1),q(:,2),'kx')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q]=cyc_shift(p,i0,type);

q=p;
n=size(p,1); m=n-i0+1;

if i0 > 1;
  if type==1;                   % forward
     q(1:m,:)=p(i0:n,:);
     q(m+1:n,:)=p(1:i0-1,:);
  else;                         % inverse
     q(i0:n,:)=p(1:m,:);
     q(1:i0-1,:)=p(m+1:n,:);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xa,ya,sa]=spline_block_arc1(x,y);

[m,n]=size(x);

dx=diff(x);dy=diff(y);s=sqrt(dx.*dx+dy.*dy);s=cumsum(s);
if m>n; s=[0;s]; else; s=[0,s]; end; 

n=length(s); nn=10*(n-1); ss=0:nn; ss=s(n)*ss/nn; ns=length(ss); 

aa=arclength(s,x,y,ss); aa=(n-1)*aa./aa(ns); 

a=0:(n-1); sa = spline(aa,ss,a); % Uniform sample in archlength;

xa=spline(s,x,sa)'; ya=spline(s,y,sa)'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xa,ya,sa]=spline_block_arc2(x,y,nout);

[m,n]=size(x);

dx=diff(x);dy=diff(y);s=sqrt(dx.*dx+dy.*dy);s=cumsum(s);
if m>n; s=[0;s]; else; s=[0,s]; end; 

n=length(s); nn=10*(n-1); ss=0:nn; ss=s(n)*ss/nn; ns=length(ss); 

aa=arclength(s,x,y,ss); aa=(nout-1)*aa./aa(ns); 

a=0:(nout-1); sa = spline(aa,ss,a); % Uniform sample in archlength;

xa=spline(s,x,sa)'; ya=spline(s,y,sa)'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [aa]=arclength(s,x,y,ss);


xx=spline(s,x,ss); yy=spline(s,y,ss);
dx=diff(xx);dy=diff(yy);aa=dx.*dx+dy.*dy;aa=sqrt(aa);
[n1,n2]=size(aa);
if n1>n2; aa=[0;aa]; else; aa=[0,aa]; end; aa=cumsum(aa);
