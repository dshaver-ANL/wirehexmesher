function d=dline(p,x1,y1,x2,y2,xc,yc)

vx = x2-x1; 
cx = xc-x1;
vy = y2-y1; 
cy = yc-y1;

vxc = vx*cy-vy*cx;

nx = -vxc*vy;
ny =  vxc*vx;

nn = sqrt(nx*nx+ny*ny);
nx = nx/nn;
ny = ny/nn;

d = nx*(p(:,1)-x1) + ny*(p(:,2)-y1);
