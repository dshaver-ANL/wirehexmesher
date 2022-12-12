
%
%   Generate only essential blocks; use 19 pin case as template
%
%


function [pin_blcko,pin_block,type_block,xyblock]=gen_ess_blocks_L(thw,kgraph,Nx,Col,Row,Rowdist,z,ne,T,S) 

hold off;


if T>0
    load wire_quad_mesh_flat; x=xx; y=yy; t=tt; s=ssl; 
else
    load wire_quad_mesh; x=xx; y=yy; t=tt; s=ssl; 
end

D = 1; P = PD*D; R=D/2;


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Generate matrix of pin coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

hold off; plot(0,0,'g.')
[npin,X,Y,xm,ym,km] = pin_coors(ne,P); 

%npin

for k=1:npin; txtstrg(X(k),Y(k),D,sprintf('%d',k));axis equal;hold on;end;
%pause(1)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Construct cells connecting pins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


[nvtx,nevtx,ntcl,tcell,necl,ecell,nccl,ccell,X,Y]=...
                       get_cells(ne,npin,km,xm,ym,X,Y,D,P);
ncell=ntcl+necl+nccl;cell=zeros(4,ncell);cell_type=zeros(ncell,1);
k0=1;k1=ntcl;k2=k1+necl;k3=k2+nccl;
cell(1:3,   1:k1)=tcell(1:3,1:ntcl); cell_type(   1:k1)=1;
cell(1:4,k1+1:k2)=ecell(1:4,1:necl); cell_type(k1+1:k2)=2;
cell(1:4,k2+1:k3)=ccell(1:4,1:nccl); cell_type(k2+1:k3)=3;

X(nvtx+1)=X(nvtx-5); Y(nvtx+1)=Y(nvtx-5);
plot(X(nvtx-5:nvtx+1),Y(nvtx-5:nvtx+1),'r-')
plot(X(nvtx-5:nvtx+1),1.0*Y(nvtx-5:nvtx+1),'r.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Define Cell centroids + additional cell vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Xcl,Ycl]=cell_center(npin,X,Y,ntcl,tcell,necl,ecell,nccl,ccell,R);
%plot(Xcl(1:ncell),Ycl(1:ncell),'ko'); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Couple cells to pins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[p2cell,ncp] = pin_2_cell(nvtx,ntcl,tcell,necl,ecell,nccl,ccell);
[cell2p,npc] = cell_2_pin(p2cell,ncp,ncell,npin);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Generate pin-based block mesh.
%
%    ve = vertex pointers for each element, e=1:nel
%    ce = curve side info for each element
%    be = boundary conditions for each element
%    pe = pin id for each element
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[nel,nelp,nele,npts,ve,ce,be,pe,qe,te,Xcl,Ycl,neigh_e,neigh_s]=...
  gen_block(ncell,npin,ne,ncp,p2cell,X,Y,Xcl,Ycl,R,ntcl,necl,ecell,ccell,T);

  hold off
  for e=1:nel;
    txt_elem(Xcl,Ycl,ve(:,e),ce(:,e),sprintf('%d',e),'b-');hold on;axis equal;axis off
   end; 
 nek_mesh_out(Xcl,Ycl,ve,nel)
 %fprintf('AFTER gen_block, B')


%
%  Shift cell centers (Xcl,Ycl) --> (Xsl, Ysl) to local maxima of fdpnwr()
%

Xsl=Xcl; Ysl=Ycl; 

fr=1.; w= [ thw , Dw , fr , R, T, S ]; P2=P/2;  % Wire definition


j=0; 
for k=1:ntcl;
  xc0=[Xcl(j+k),Ycl(j+k)];
  xcm=find_max_d(tcell(:,k),Xcl(j+k),Ycl(j+k),X,Y,R,P,w,npin,nevtx,xc0); 
  Xsl(j+k)=xcm(1); Ysl(j+k)=xcm(2);plot(xcm(1),xcm(2),'m*');
  %pause
end;


j=ntcl; 
for k=1:necl;
  xc0=[Xcl(j+k),Ycl(j+k)];
  xcm=find_max_d(ecell(:,k),Xcl(j+k),Ycl(j+k),X,Y,R,P,w,npin,nevtx,xc0); 
  Xsl(j+k)=xcm(1); Ysl(j+k)=xcm(2);plot(xcm(1),xcm(2),'ro');
end;

j=ntcl+necl; 
for k=1:nccl; cl = j+k;
  xc0=[Xcl(j+k),Ycl(j+k)]'; nn=xc0/norm(xc0);

  xcm=xy_proj(xc0,nn,cl,cl,X,Y,cell,cell_type,cell2p,npc,P,w,nevtx);
% xcm=find_max_d(ccell(:,k),Xcl(j+k),Ycl(j+k),X,Y,R,P,w,npin,nevtx,xc0); 
  Xsl(j+k)=xcm(1); Ysl(j+k)=xcm(2);plot(xcm(1),xcm(2),'ro');
end;
  %pause(.01);

  

  
[nel,nelp,nele,npts,ve,ce,be,pe,qe,te,Xsl,Ysl,neigh_e,neigh_s]=...
  gen_block(ncell,npin,ne,ncp,p2cell,X,Y,Xsl,Ysl,R,ntcl,necl,ecell,ccell);


%%for e=1:nel;
%%  txt_elem(Xsl,Ysl,ve(:,e),ce(:,e),sprintf('%d',e),'k-');hold on;axis equal;axis off
%%end; 


%epin
%  Connect shifted cell centers
%



flag = spalloc(4*nel,4*nel,4*nel); nedge=0;
for e=1:nel;
   c3=ve(3,e);c4=ve(4,e);

    id=pe(e);
    if id>0
        thc3=atan2(Ysl(c3)-Y(id),Xsl(c3)-X(id));
        thc4=atan2(Ysl(c4)-Y(id),Xsl(c4)-X(id));
    else
        thc3=atan2(Ysl(c3),Xsl(c3));
        thc4=atan2(Ysl(c4),Xsl(c4));
    end
    
    if (thc3>thc4 && abs(thc3-thc4)<pi)||(thc3<thc4 && abs(thc3-thc4)>pi)
        C1=c4; C2=c3;
    else
        C1=c3; C2=c4;
    end
   
   if flag(C1,C2)==0; nedge=nedge+1;flag(C1,C2)=nedge;
      %plot([Xsl(c34),Xsl(C34)],[Ysl(c34),Ysl(C34)],'r-')
   end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


medge=Col*Nx+1;
m2=Row*Nx+1;

% Distributions for sides and top
distC=zeros(Col*Nx+1,1);
n=0;
for i=1:Col 
    for j=1:Nx
        n=n+1;
        distC(n)=(i-1)/Col + z(j)/Col;
    end
end
distC(n+1)=1;

distR=zeros(Row*Nx+1,1);
n=0;
r=0;
for i=1:Row
   for j=1:Nx
      n=n+1;
      distR(n)=r+Rowdist(i)*z(j)/Row;
   end
   r=r+Rowdist(i)*z(Nx+1)/Row;
end
distR(n+1)=1;




m=medge; xedge=zeros(m,4*nel); yedge=xedge; dedge=xedge; flag=0.*flag; nedge=0;
for e=1:nel; 
    
    c3=ve(3,e);c4=ve(4,e);

% Pin section points go ccw around the pins, wall section points go
% ccw around the center
    id=pe(e);
    if id>0
        thc3=atan2(Ysl(c3)-Y(id),Xsl(c3)-X(id));
        thc4=atan2(Ysl(c4)-Y(id),Xsl(c4)-X(id));
    else
        thc3=atan2(Ysl(c3),Xsl(c3));
        thc4=atan2(Ysl(c4),Xsl(c4));
    end
    
    if (thc3>thc4 && abs(thc3-thc4)<pi)||(thc3<thc4 && abs(thc3-thc4)>pi)
        C1=c4; C2=c3;
    else
        C1=c3; C2=c4;
    end
    

    if flag(C1,C2)==0; nedge=nedge+1;flag(C1,C2)=nedge;
        
        pc=qe(e);
        
      [xedge(:,nedge),yedge(:,nedge),dedge(:,nedge)]=...
            xy_edge(C1,C2,Xsl,Ysl,X,Y,cell,cell_type,cell2p,npc,P,w,m,nevtx,distC,pc,xedge,yedge,T,S);
        
        %plot(xedge,yedge,'b*')
        %pause
        %{
        plot(Xsl(C1),Ysl(C1),'m*')
 
        plot(Xsl(C2),Ysl(C2),'r*')
        plot(X(pc),Y(pc),'kx')
        

        pause
        plot(X(pc),Y(pc),'yx')
        plot(Xsl(C1),Ysl(C1),'y*')
        plot(Xsl(C2),Ysl(C2),'y*')
        
        
        pause
   %}     
    end;
    
end;

%pause

%size(Xsl)
%size(Ysl)
%size(X)
%size(Y)

%plot(Xsl,Ysl,'bx')
%plot(X,Y,'r*')

%size(xedge)
%pause
%plot(xedge,yedge,'gx')
%pause

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Draw pin + wire boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%



pin=1; wt= w; [ pi/2 , Dw , fr , R, T ]; % Wire definition
th0 = 0;X0=[X(pin)+R*cos(th0),Y(pin)+R*sin(th0)];
mp = 200; h = 1.05569933892436*(2*pi*R/mp); th=1; iter=0;
x0=X0; xpin=zeros(mp+1,2); xpin(1,:)=x0;
for k=1:mp;
  [x2,x1]=pin_tan_step(x0,h,X,Y,pin,wt); xpin(k+1,:)=x2; x0=x2;
end;

%plot(xpin(:,1),xpin(:,2),'rx')
%pause

xpn=xpin;
for pin=2:npin;
   xpn(:,1)=xpin(:,1) + X(pin)-X(1);
   xpn(:,2)=xpin(:,2) + Y(pin)-Y(1);
  %plot(xpn(:,1),xpn(:,2),'b-'); 
end;
%pause
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    BLOCK GENERATION:  For each block, connect to pin or edge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Cells are determined by c3 and c4.   Pin is determined by pe(e)
%















% adds nodes to outside walls

ct_max = 0;
for e=1:nel; c1=ve(1,e);c4=ve(4,e);pin=qe(e);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  if flag(c1,c4)==0; nedge=nedge+1;flag(c1,c4)=nedge;
    if pin>0;
      [xedge(1:m,nedge),yedge(1:m,nedge)]=xy_pin(c1,c4,e,pin,Xsl,Ysl,X,Y,w,m,T);
    else  % Corner or Edge cell
       c2=ve(2,e); c3=ve(3,e);
       ct_max        = max(ct_max,cell_type(c4));
       cell_type(c4) = ct_max;

       m1=m-1; s = (0:m1)/m1;
       xedge(1:m,nedge) = Xsl(c4) + s.*(Xsl(c1)-Xsl(c4));
       yedge(1:m,nedge) = Ysl(c4) + s.*(Ysl(c1)-Ysl(c4));
    end;
  end;
  %plot(xedge(1:m,nedge),yedge(1:m,nedge),'bx'); 
  %plot(xedge,yedge,'bx');
  %pause
end;



%%%  Count active elements and identify pin_block and typ_block

%epin = [ 1 3 8 13 11 6 2 5  10 12 9 4];
%cpin = [ 1 2 3 4 0 0 5 6 0 7 0 8 9 0 0 10 11 12 13];
[cpin, epin]=make_arrays(ne);

pin_block=zeros(1,30*ne-24);
pin_blcko=pin_block;
type_block=pin_block;


ee = 0;  
for e=1:nel;
  pin=qe(e);
  if pin>0 && cpin(pin)>0
     ee             = ee+1;
     pin_block (ee) = cpin(pin);
     pin_blcko (ee) = pin;
     type_block(ee) = te(e);

  elseif pin < 0;

     ee             = ee+1;
     pin_block(ee)  = epin(abs(pin));
     pin_blcko(ee)  = pin;
     type_block(ee) = te(e);

  end;
end;



nel_ess = ee; ee2e = zeros(nel_ess,1); e2ee=zeros(nel,1);
 

xyblock=zeros(m2,medge,nel_ess,2); 

%tic;

if T>0
    load wire_quad_mesh_flat; x=xx; y=yy; t=tt; s=ssl; 
else
    load wire_quad_mesh; x=xx; y=yy; t=tt; s=ssl; 
end

ee = 0; pinlast = 0;
for e=1:nel;    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  c3 =ve(3,e); c4 =ve(4,e); 
  
  % Find lower angle of the two
  id=pe(e);
    if id>0
        thc3=atan2(Ysl(c3)-Y(id),Xsl(c3)-X(id));
        thc4=atan2(Ysl(c4)-Y(id),Xsl(c4)-X(id));
    else
        thc3=atan2(Ysl(c3),Xsl(c3));
        thc4=atan2(Ysl(c4),Xsl(c4));
    end
    
    if (thc3>thc4 && abs(thc3-thc4)<pi)||(thc3<thc4 && abs(thc3-thc4)>pi)
        C1=c4; C2=c3;
    else
        C1=c3; C2=c4;
    end
  
  c12=ve(2,e);C12=ve(1,e); if C1==c4; c12=ve(1,e);C12=ve(2,e); end;

  kedge = flag(C1,C2); pin=qe(e); 

  e2ee(e) = 0;

   if pin>0 && cpin(pin)>0 
    ee=ee+1; ee2e(ee) = e; e2ee(e) = ee;
   
      tp=0;
      change=0;
      for i=1:medge; 
         tpa=tp;
         [xx,yy,tp]=xy_pin2(pin,xedge(i,kedge),yedge(i,kedge),X,Y,w,m2,distR,x,y,t,s);
         
         xyblock(:,i,ee,1)=flipud(xx); xyblock(:,i,ee,2)=flipud(yy);
         if tpa~=tp && i>1; change=i;
         end;
      end;
      
      
      %plot(xedge,yedge,'rx')
      
      %plot(xyblock(1,:,ee,1),xyblock(1,:,ee,2),'r*')
      %plot(xyblock(end,:,ee,1),xyblock(end,:,ee,2),'k*')
      %plot(xyblock(:,:,ee,1),xyblock(:,:,ee,2),'bx')
      
      %if ee>1
       %   plot(xyblock(1,:,ee-1,1),xyblock(1,:,ee-1,2),'bx')
      %end
      
      %pause
      
      %%%%%%%%%%%%%% Landon
      %Check if this element has bad radial lines at a squeezed portion
      %{           
      NOt needed
      
      if tp==2  % only for type two blocks
         in=squeeze(xyblock(1,:,ee,:))
         out=squeeze(xyblock(end,:,ee,:))
         size(in)
         size(out)
         
         for iout=1:medge
            [mdist(iout),iin(iout)]=min(((out(iout,1)-in(:,1)).^2+(out(iout,2)-in(:,2)).^2).^(1/2))
    
            plot(out(iout,1)+(0:0.1:1)*(in(iin(iout),1)-out(iout,1)),out(iout,2)+(0:0.1:1)*(in(iin(iout),2)-out(iout,2)),'k')
         end

         [md,miout]=min(mdist)
         plot(out(miout,1)+(0:0.1:1)*(in(iin(miout),1)-out(miout,1)),out(miout,2)+(0:0.1:1)*(in(iin(miout),2)-out(miout,2)),'m')
         
          
      end
      pause
      %}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
  
     
      % Fix spacing of certain elements
      if change>0
         if change==medge; change=change-1; end;
         conelem=ceil((change-1)/(Nx));
         conindex=(conelem-1)*(Nx)+(1:Nx+1);
         th=atan2(xyblock(1,conindex,ee,2)-Y(pin_blcko(ee)),xyblock(1,conindex,ee,1)-X(pin_blcko(ee)));
         rad=sqrt((xyblock(1,conindex,ee,1)-X(pin_blcko(ee))).^2+(xyblock(1,conindex,ee,2)-Y(pin_blcko(ee))).^2);
         xyblock(1,conindex,ee,1)=rad.*cos(th(1)+z'*(th(Nx+1)-th(1)))+X(pin_blcko(ee));
         xyblock(1,conindex,ee,2)=rad.*sin(th(1)+z'*(th(Nx+1)-th(1)))+Y(pin_blcko(ee));
         
       
         
         for i=conindex
            xyblock(:,i,ee,1)=xyblock(1,i,ee,1)+flipud(1-distR).*(xyblock(end,i,ee,1)-xyblock(1,i,ee,1));
            xyblock(:,i,ee,2)=xyblock(1,i,ee,2)+flipud(1-distR).*(xyblock(end,i,ee,2)-xyblock(1,i,ee,2));
         end
         
         
         
        
         %if conindex(1)==1
         %   fprintf('done \n')
         %   xyblock(:,1,ee,1)=(xyblock(:,end,ee-1,1)+xyblock(:,1,ee,1))/2;
         %   xyblock(:,1,ee,2)=(xyblock(:,end,ee-1,2)+xyblock(:,1,ee,2))/2;
         %   xyblock(:,end,ee-1,1)=xyblock(:,1,ee,1);
         %   xyblock(:,end,ee-1,2)=xyblock(:,1,ee,2);
         %end
         
         
         
  
      end

      if ee>1 && pin_blcko(ee)==pin_blcko(ee-1)
         xyblock(:,1,ee,1)=(xyblock(:,end,ee-1,1)+xyblock(:,1,ee,1))/2;
         xyblock(:,1,ee,2)=(xyblock(:,end,ee-1,2)+xyblock(:,1,ee,2))/2;
         xyblock(:,end,ee-1,1)=xyblock(:,1,ee,1);
         xyblock(:,end,ee-1,2)=xyblock(:,1,ee,2);
      end

      efirst=ee-nnz(pin_blcko==pin_blcko(ee))+1;
      if ee>4 && pin_blcko(ee)==pin_blcko(efirst)
         xyblock(:,end,ee,1)=(xyblock(:,1,efirst,1)+xyblock(:,end,ee,1))/2;
         xyblock(:,end,ee,2)=(xyblock(:,1,efirst,2)+xyblock(:,end,ee,2))/2;
         xyblock(:,1,efirst,1)=xyblock(:,end,ee,1);
         xyblock(:,1,efirst,2)=xyblock(:,end,ee,2);    
      end
      
      %{
      if ee>1
        plot(xyblock(:,:,ee-1,1),xyblock(:,:,ee-1,2),'mx')
        plot(xyblock(:,end,ee-1,1),xyblock(:,end,ee-1,2),'m*')
      end
      plot(xyblock(:,:,ee,1),xyblock(:,:,ee,2),'bx')
      plot(xyblock(:,1,ee,1),xyblock(:,1,ee,2),'b*')
      pause
      %}

    

    if mod(kgraph,4)==0; 
        %spline_block(xyblock(:,:,ee,1),xyblock(:,:,ee,2)); 
    end;

  elseif pin < 0; ee=ee+1; ee2e(ee) = e; e2ee(e) = ee;


    x12=Xsl(c12); X12=Xsl(C12); y12=Ysl(c12); Y12=Ysl(C12);
%   plot(x12,y12,'go');plot(X12,Y12,'y*');
    %m21=m2-1; s = (0:m21)/m21;
    s=distR';

    [amin,imin]=min(dedge(:,kedge));

    if imin==1 || imin==medge; jj=0;   % If min pt is at start or end.
       
       for i=1:medge; jj=jj+1; 
         xtarg = x12+(i-1)*(X12-x12)/(medge-1);
         ytarg = y12+(i-1)*(Y12-y12)/(medge-1);
         xyblock(:,jj,ee,1)=xedge(i,kedge)+s.*(xtarg-xedge(i,kedge));
         xyblock(:,jj,ee,2)=yedge(i,kedge)+s.*(ytarg-yedge(i,kedge));
       end;
    else  % Min pt is in middle of element
      
       xt=zeros(medge);
       yt=xt;

       xt=X12*distC+x12*(1-distC);
       yt=Y12*distC+y12*(1-distC);
       
       
       for i=1:medge; 
         xtarg = xt(i); ytarg = yt(i);
         xyblock(:,i,ee,1)=xedge(i,kedge)+s.*(xtarg-xedge(i,kedge));
         xyblock(:,i,ee,2)=yedge(i,kedge)+s.*(ytarg-yedge(i,kedge));
       end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Check if this element has bad radial lines at a squeezed portion
      
      in=squeeze(xyblock(1,:,ee,:));
      out=squeeze(xyblock(end,:,ee,:));

      %plot(in(:,1),in(:,2),'rx')
      %plot(out(:,1),out(:,2),'bx')
        
      for iin=1:medge
         [mdist(iin),iout(iin)]=min(((in(iin,1)-out(:,1)).^2+(in(iin,2)-out(:,2)).^2).^(1/2));
    
         %plot(in(iin,1)+(0:0.1:1)*(out(iout(iin),1)-in(iin,1)),in(iin,2)+(0:0.1:1)*(out(iout(iin),2)-in(iin,2)),'k')
         
      end
      

      [md,miin]=min(mdist);    
      
      xin=[in(miin,1),in(miin,2)];
      xout=[out(iout(miin),1),out(iout(miin),2)];
      
      
      if miin>1 && miin<medge
      
      
      low=out(iout(miin)-1,:);
      high=out(iout(miin)+1,:);
      
      if iout(miin)==1; low=out(iout(miin),:); end
      if iout(miin)==medge; high=out(iout(miin),:); end
      
      if high(1)==low(1); high(1)=low(1)+1E-6; end
      if high(2)==low(2); high(2)=low(2)+1E-6; end
      
      
      %{
      xout(1)
      xout(1)-(low(1):(high(1)-low(1))/100:high(1))
      (xout(1)-(low(1):(high(1)-low(1))/100:high(1))).^2
      xout(2)
      low(2)
      high(2)
      (xout(2)-(low(2):(high(2)-low(2))/100:high(2))).^2
      ( (xout(1)-(low(1):(high(1)-low(1))/100:high(1))).^2+(xout(2)-(low(2):(high(2)-low(2))/100:high(2))).^2 ).^(1/2)
      %}
      
      
      
      [mdr,index]=min( ( (xin(1)-(low(1):(high(1)-low(1))/100:high(1))).^2+(xin(2)-(low(2):(high(2)-low(2))/100:high(2))).^2 ).^(1/2) );
      
      %pause
      
      xforeal=low+index/100*(high-low);
      
      xout=xforeal;
      
      
          
      %{
      if iout(miin)>1
          [mcw,icw]=min(((xin(1)-(1:-0.1:0)*xout(1)+(0:0.1:1)*out(iout(miin)-1,1)).^2+(xin(2)-(1:-0.1:0)*xout(2)+(0:0.1:1)*out(iout(miin)-1,2)).^2).^(1/2))
          xmcw=icw/10*xout+(1-icw/10)*out(iout(miin)-1,:);
  
      %}
      
      %xout=[out(iout(miin),1),out(iout(miin),2)];
      %plot(in(miin,1)+(0:0.1:1)*(out(iout(miin),1)-in(miin,1)),in(miin,2)+(0:0.1:1)*(out(iout(miin),2)-in(miin,2)),'m')
      
      %plot(xin(1),xin(2),'b*')
      %plot(xout(1),xout(2),'r*')
      
      %
      % if the length of line from in(miin) to out(miin) is 2X longer than
      % md then do this:
      %dist=((xin(1)-out(miin,1))+(xin(2)-out(miin,2)))^(1/2);
      %if dist > 2*md
         
          xt(1:miin)=out(1,1)+(0:1/(miin-1):1)*(xout(1)-out(1,1));
          yt(1:miin)=out(1,2)+(0:1/(miin-1):1)*(xout(2)-out(1,2));  
          
          xt(miin:medge)=xout(1)+(0:1/(medge-miin):1)*(out(medge,1)-xout(1));
          yt(miin:medge)=xout(2)+(0:1/(medge-miin):1)*(out(medge,2)-xout(2));

          
          for i=1:medge; 
             xtarg = xt(i); ytarg = yt(i);
             xyblock(:,i,ee,1)=xedge(i,kedge)+s.*(xtarg-xedge(i,kedge));
             xyblock(:,i,ee,2)=yedge(i,kedge)+s.*(ytarg-yedge(i,kedge));
          end;
      %end
          
      %plot(xyblock(:,:,ee,1),xyblock(:,:,ee,2),'bx')
      %plot(xout(1),xout(2),'r*')
      %plot(high(1),high(2),'g*')
      %plot(low(1),low(2),'k*')
      %pause
      end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end;
%   toc; 'PIN < 0'
   
     

    %if mod(kgraph,4)==0; 
    spline_block(xyblock(:,:,ee,1),xyblock(:,:,ee,2)); 
     %end;
  end;
  if pin~=pinlast; drawnow; pause(.01); end; pinlast = pin;
end;

if kgraph==0;
  save data_ess_cell ...
       Xcl Ycl nel nelp nele npts ve ce be pe qe te neigh_e neigh_s ...
       X Y ntcl tcell necl ecell nccl ccell R npin ne ncp p2cell e2ee ee2e
end;




%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    END OF gen_ess_blocks.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Analysis Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%{
function [o,rad]=rad_circ3(xc) % Given 3 points, find radius of circle %1


plot(xc(:,1),xc(:,2),'ko')
x1=xc(1,:)'; x2=xc(2,:)'; x3=xc(3,:)'; x21=x2-x1; x32=x3-x2; x31=x3-x1;
t1=x21/norm(x21); t2=x32/norm(x32); n1=[-t1(2),t1(1)]'; n2=[-t2(2),t2(1)]';
A = [ n1 , -n2 ]; Ai=inv(A); alpha = .5*Ai*x31;
o = alpha(1)*n1 + .5*(x1+x2); r = o-x2; rad = norm(r);

%}
%txtcirc(o(1),o(2),2*rad,'curve')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=get_curve(x1,x2,xcg,rad,nc) %2


xm=.5*(x1+x2); l12=norm(x2-x1); th=(x2-x1)/l12; nh=[ -th(2) , th(1) ]';

dcg = xcg-xm; if nh'*dcg < 0; nh=-nh; end;  %  Make vector point to center

if rad < 0; nh=-nh; end;                    %  Flip if concave element

h = rad*rad - .25*l12*l12; h=sqrt(h);
O = xm + h*nh;

dt = -.5*l12/abs(rad); dt=2*asin(dt)/(nc+1); if rad>0; dt=-dt; end;
th0=atan2(x1(2)-O(2),x1(1)-O(1)); th=th0 + dt*(1:nc)'; 
x = O(1) + abs(rad)*cos(th); y = O(2) + abs(rad)*sin(th);



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Plotting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function txtstrg(x,y,D,txt) %3

text(x-.10,y,txt);
axis equal; axis off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function txtcirc(x,y,D,txt) %4

R=D/2; th=pi*(0:.05:2);
plot(x+R*cos(th),y+R*sin(th),'r-'); text(x-.10,y,txt);
axis equal; axis off
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt_elem(x,y,v,c,txt,SS) %5


cmx=max(c); cmn=min(c);  % Any curve sides

if cmx==0 && cmn==0;      % No curve sides
   v4 = [ v ; v(1) ]; plot(x(v4),y(v4),'b-');
   xb = sum(x(v))/length(v); yb = sum(y(v))/length(v); text(xb-.10,yb,txt);
   axis equal; axis off
else
   nc=20; xx=zeros(4*nc+2,1);yy=xx; xcg = [sum(x(v))/4, sum(y(v))/4]';
   m=0; 
   for k=1:4;
      m=m+1; xx(m)=x(v(k)); yy(m)=y(v(k));
      if abs(c(k))>0;
	 k1 = k+1; if k1>4; k1=1; end;
         x1 = [x(v(k)),y(v(k))]'; x2 = [x(v(k1)),y(v(k1))]';
         [xx(m+1:m+nc),yy(m+1:m+nc)]=get_curve(x1,x2,xcg,c(k),nc); m=m+nc;
      end;
   end;
   m=m+1; xx(m)=x(v(1)); yy(m)=y(v(1));
   plot(xx(1:m),yy(1:m),SS); 
   xb = sum(x(v))/length(v); yb = sum(y(v))/length(v); text(xb-.10,yb,txt);
end;
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function txt_cell(x,y,v,txt) %6


v4 = [ v ; v(1) ]; plot(x(v4),y(v4),'b-');
xb = sum(x(v))/length(v); yb = sum(y(v))/length(v); text(xb-.10,yb,txt);
axis equal; axis off
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nek_mesh_out(x,y,cell,nel) %7


fid=fopen('mesh.1','w');
fprintf(fid,'%8i   %8i   %8ig\n',nel,2,nel);
for e=1:nel;
   fprintf(fid,'            ELEMENT%5i [    1A]    GROUP     0\n',e);
   fprintf(fid,'%15.7g  %15.7g  %15.7g  %15.7g\n',x(cell(:,e))');
   fprintf(fid,'%15.7g  %15.7g  %15.7g  %15.7g\n',y(cell(:,e))');
end;
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdpins(p,xc,yc,rad) %8


d = dcircle(p,xc(1),yc(1),rad); nc=length(xc);
for i=2:nc;
   g = dcircle(p,xc(i),yc(i),rad); 
   d = min(d,g);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdedge(p,xe,ye) %9

d=0*p(:,1) + 1.e20; ne=length(xe); 

xc=sum(xe)/ne; yc=sum(ye)/ne;
x1=xe(ne);     y1=ye(ne);


for i=1:ne;
   x0=x1;y0=y1; x1=xe(i); y1=ye(i);
   g = dline(p,x0,y0,x1,y1,xc,yc); 
   d=min(d,g);
end;


%{
xe2=xe;ye2=ye;
for i=1:ne
   xe2(i)=x1;ye2(i)=y1; x1=xe(i); y1=ye(i);
end
g=dline2(p,xe2,ye2,xe,ye,xc,yc);
d=min(g,[],2);
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grd=pin_grad(x,xc,yc,wire) %10


h=1.e-5;
px = [ x(1)-h , x(1)+h ,  x(1)   , x(1)   ]'; 
py = [ x(2)   , x(2)   ,  x(2)-h , x(2)+h ]'; 
X = [ px , py ]; Z = fdpnwr(X,xc,yc,wire);
grd=0*x; grd(1)=.5*(Z(2)-Z(1))/h; grd(2)=.5*(Z(4)-Z(3))/h;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdpnwr(p,xc,yc,wire) %11


  R = wire(4);  % Pin radius
  d=fdpins(p,xc,yc,R); dw=fdwire(p,xc,yc,wire,0); d=min(d,dw);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Modified Wire
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function d=fdwire_big_circle(p,xc,yc,wire) %12


th  = wire(1);  % Angle of wire
dw  = wire(2);  % wire diameter
rf  = wire(3);  % fillet ratio
rad = wire(4);  % Pin radius

D   = 2*rad;   DD = D+dw;  RR=DD/2;
rw  = dw/2;

xw=xc+rw*cos(th); yw=yc+rw*sin(th);

nc=length(xc)
for k=1:nc;
  if k==1; d = dcircle(p,xw(k),yw(k),RR); else;
     g = dcircle(p,xw(k),yw(k),RR); d = min(d,g);
  end;
end;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Standard Wire
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdwire(p,xc,yc,wire,tf) %13









th  = wire(1);  % Angle of wire
dw  = wire(2);  % wire diameter
rf  = wire(3);  % fillet ratio
rad = wire(4);  % Pin radius
T = wire(5);
S = wire(6);

rw = .5*dw; rf=rf*rw;

xw=xc+(rad+rw-S)*cos(th); yw=yc+(rad+rw-S)*sin(th);

thf = ( (rad+rw-S).^2+(rad+rf).^2-(rw+rf).^2 ) / (2.*(rad+rf)*(rad+rw-S));
thf = acos(thf);

xtri=zeros(3,1); ytri=zeros(3,1); 

nc=length(xc);
for k=1:nc;
  if k==1; d = dcircle(p,xw(k),yw(k),rw); else
     g = dcircle(p,xw(k),yw(k),rw); d = min(d,g);
  end;
  
  

    %%%%%%%%%%%%%%%%%%%%%%%% Modified here Landon
  
  if (T>0)
      
     
    line_x=rad+2*rw-T-S;
    xp=p(:,1)-xc(k);
    yp=p(:,2)-yc(k);
    dist=(xp.^2+yp.^2).^(1/2);
    th2=atan2(yp,xp);
    th2=th2-th;
  
    %these=(p(:,1)>(line_x-T) & abs(p(:,2))<(rw-T)*tan(acos((rw-T)/rw)));
    these=(dist.*cos(th2)>(line_x-2*T) & abs(dist.*sin(th2))<(rw-T)*tan(acos((rw-T)/rw)));
    
    di=dist.*cos(th2)-line_x;
    d(these)=di(these);

  end

  
  
  
  
  
  
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %}

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
%pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xcm=find_max_d(cell,Xcl,Ycl,X,Y,R,P,wire,npin,nevtx,xcm0) %14


ipin  = find(cell < npin+1);
cellp = cell(ipin);
th    = wire(1);  % Angle of wire
dw    = wire(2);  % wire diameter
rf    = wire(3);  % fillet ratio
R     = wire(4);  % Pin radius
P2    = P/2    ;  % Pitch / 2

h = dw/3; xcm=xcm0; maxiter=20; dmxx=0; xmxx=xcm;
for k=1:maxiter;
   nx=4; if k==1; nx=9; end;
   h=3*h/4;
   xmn=xcm(1)-h; xmx=xcm(1)+h; ymn=xcm(2)-h; ymx=xcm(2)+h;
   dx=h/nx; dy=dx; xx=xmn:dx:xmx; yy=ymn:dy:ymx;
   [XX,YY]=meshgrid(xx,yy); nx=size(XX,1); ny=size(XX,2);
   Xt=reshape(XX,nx*ny,1); Yt=reshape(YY,nx*ny,1); PP=[ Xt , Yt ];

   vn = [ cos(th) , sin(th) ];        % Starting at cell center, search only in 
   vp = PP-ones(nx*ny,1)*[Xcl,Ycl];   % wedge that is aligned with pin angle
   iv = find( vp*vn' > 0);
   PP=PP(iv,:);

   dp=fdpins(PP,X(cellp),Y(cellp),R);
   dw=fdwire(PP,X(cellp),Y(cellp),wire,0);
   de=fdedge(PP,X(nevtx+1:nevtx+6),Y(nevtx+1:nevtx+6));
   d=min(dp,de); d=min(d,dw); d=max(d,0.); 
   dp = max(dp,0); de=max(de,0); dw=max(dw,0); dp=dp.*dw.*de;
   d=d - .001*dp;
   [dmx,imx]=max(d);
   xcm(1)=PP(imx,1); xcm(2)=PP(imx,2);
   if dmx > dmxx; dmxx=dmx; xmxx=xcm; end; % Keep max point
end;
xcm=xmxx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xedge,yedge,dm]=...
            xy_edge(c1,c2,Xcl,Ycl,X,Y,cell,cell_type,cell2p,npc,P,wire,m,nevtx,dist,pc,xedgec,yedgec,T,S) %15
        

plot(X(pc),Y(pc),'rx')

 th    = wire(1);  % Angle of wire
 dw    = wire(2);  % wire diameter
% rf    = wire(3);  % fillet ratio
 R     = wire(4);  % Pin radius
% P2    = P/2    ;  % Pitch / 2

%
%  Now, march along the line connecting c1 to c2 and find projection
%  onto locus of points equidistant from cell boundaries
%

X1 = [Xcl(c1),Ycl(c1)]'; X2 = [Xcl(c2),Ycl(c2)]';
t=(X2-X1);tn=t/norm(t);nn=[-tn(2),tn(1)]';

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%% Landon

%plot(X1(1),X1(2),'g*')
%plot(X2(1),X2(2),'r*')

eps=0.001;

%p1b=intersect(find(abs(xedgec(1,:)-X1(1))<eps),find(abs(yedgec(1,:)-X1(2))<eps))
p1e=intersect(find(abs(xedgec(end,:)-X1(1))<eps),find(abs(yedgec(end,:)-X1(2))<eps));
p2b=intersect(find(abs(xedgec(1,:)-X2(1))<eps),find(abs(yedgec(1,:)-X2(2))<eps));
%p2e=intersect(find(abs(xedgec(end,:)-X2(1))<eps),find(abs(yedgec(end,:)-X2(2))<eps))
connected=0;
int=intersect(p1e,p2b);


Xc=X(pc); Yc=Y(pc);
X0=zeros(2,m);
X0(:,1)=X1;
X0(:,m)=X2;
dX=norm(X2-X1)/(m-1);


dx1=X1(1)-Xc; dy1=X1(2)-Yc;
dx2=X2(1)-Xc; dy2=X2(2)-Yc;
dx12=X2(1)-X1(1); dy12=X2(2)-X1(2);
d1=(dx1^2+dy1^2)^(1/2);
%d2=((X2(1)-Xc)^2+(X2(2)-Yc)^2)^(1/2);
d12=(dx12^2+dy12^2)^(1/2);
th1=atan2(dy1,dx1);
th2=atan2(dy2,dx2);
th12=atan2(dy12,dx12);
%pi or 180???
lsf=1.5;
if (th>=th1 && th<=th2)
    tha=th-th1;
    thb=pi-th12+th1;
    thc=pi-tha-thb;
    
    di=sin(tha)*d1/sin(thc);
    xi=X1(1)+di/d12*(X2(1)-X1(1));
    yi=X1(2)+di/d12*(X2(2)-X1(2));
    Xi=[xi; yi];
    %plot(Xi(1),Xi(2),'m*')
    %pause

    [xth,dmth]=...
    xy_proj(Xi,nn,c1,c2,X,Y,cell,cell_type,cell2p,npc,P,wire,nevtx);

    gap=((xth(1)-Xc)^2+(xth(2)-Yc)^2)^(1/2)-(R+dw-T-S);
    gap=gap/(T/2);
    %if gap is smallest, smallest length is 1/2 dx, if gap is 5 times this, smallest length is 1 dx
    %lsf=(1/8*gap+3/8);
    lsf=(1/18*gap+8/18);
    mf=2-lsf;
end



if (th>=th1 && th<=th2 && lsf<1) %|| (thn>th1 && thn<th2)
ls=dX*lsf;

%plot(Xc,Yc,'k*')
%plot(xth(1),xth(2),'m*')


    %ls=dX/2;
    %mf=1.5;

    dbef=di/d12;
    daft=1-di/d12;
    nbef=round((m-1)*dbef^(1/2)/(dbef^(1/2)+daft^(1/2)));
    naft=round((m-1)*daft^(1/2)/(dbef^(1/2)+daft^(1/2)));
    nth=nbef+1;
    
    dbef=dbef*d12;
    daft=daft*d12;
  
    
    if (nbef+naft>(m-1))
        'sum of points is wrong'
        pause
    end

    eps=0.00001;
    
    % before
    ra=1.05;
    rb=1.1;
    %ls=dX/2;
    Sna = ls*(1-ra^nbef)/(1-ra)-dbef;

    
    capb=1;
    true=0;
   iter=0;
    
   %dX
   %dbef
   %tn
   %dX
   dxbef=dbef/nbef;
   
    while true==0 && nbef>0
        iter=iter+1;
            err=1;
       while abs(err)>eps
       
          Snb = ls*(1-rb^nbef)/(1-rb)-dbef;
          rc = rb-(rb-ra)*Snb/(Snb-Sna);
          Snc = ls*(1-rc^nbef)/(1-rc)-dbef;
          err=Snc;
          ra=rb;
          rb=rc;
          Sna=Snb;       
       end
       rbef=rc;
       large=ls*(rbef^nbef)/dX;
       

       if large<mf*dxbef/dX
            true=1;
       elseif iter>nbef-1
           'something aint right'
           pause
       else
           capb=capb+1;
           X0(:,capb)=X0(:,capb-1)+mf*dxbef*tn;
           nbef=nbef-1;
           dbef=dbef-mf*dxbef;
       end
    end

    %largestbefore=ls*rbef^nbef;
    %before=largestbefore/dX;
    
    
    % after
    err=1;
    ra=1.05;
    rb=1.1;
    %ls=dX/2;
    Sna = ls*(1-ra^naft)/(1-ra)-daft;
       
    capa=m;
    true=0;
    iter=0;
    
    dxaft=daft/naft;

    
    while true==0 && naft>0
       iter=iter+1;
       err=1;
       while abs(err)>eps
       
           Snb = ls*(1-rb^naft)/(1-rb)-daft;
           rc = rb-(rb-ra)*Snb/(Snb-Sna);
           Snc = ls*(1-rc^naft)/(1-rc)-daft;
           err=Snc;
           ra=rb;
           rb=rc;
           Sna=Snb;       
       end
       raft=rc;
       large=ls*(raft^naft)/dX;

%
       if large<mf*dxaft/dX
            true=1;
       elseif iter>naft-1
           'something aint right 2'
           pause
       else
           capa=capa-1;
           X0(:,capa)=X0(:,capa+1)-mf*dxaft*tn;
           naft=naft-1;
           daft=daft-mf*dxaft;
       end
    end       
%
    %largestafter=ls*raft^naft;
    %after=largestafter/dX;
   
    
    X0(:,1)=X1;
    X0(:,m)=X2;
    X0(:,nth)=[xi,yi]; 
    
    for i=nth-1:-1:capb+1
       n=-i+nth-1;
       X0(:,i)=X0(:,i+1)-ls*rbef^n*tn;
    end
        
    for i=nth+1:capa-1
        n=i-(nth+1);
        X0(:,i)=X0(:,i-1)+ls*raft^n*tn;
    end
    
    %plot(X0(1,:),X0(2,:),'bx')
    %pause

        %plot(xi,yi,'m*')
        %plot(X1(1),X1(2),'m*')
        %plot(X2(1),X2(2),'m*')
        %plot(X0(1,nbef+2:m-1),X0(2,nbef+2:m-1),'r*')
    %    plot(X1(1),X1(2)+dbef,'r*')
    %xdf=X0(1,nbef+1)+ls*tn(1)
    %ydf=X0(2,nbef+1)+ls*tn(2)
   
    %plot(xdf,ydf,'g*')
  
elseif int>0
    %connected=1;
    X0(1,:)=flipud(xedgec(:,int));
    X0(2,:)=flipud(yedgec(:,int));
    X0(:,1)=X1;
    X0(:,end)=X2;
    
    %for y=1:size(X0,2)
    %    plot(X0(1,y),X0(2,y),'m*');
    %    pause
    %nd
    %'it is so!'
    %pause
    
elseif (th1-th)<=pi/8 && (th1-th)>=0 && lsf<1
   %'below'
    err=1;
    ra=1.05;
    rb=1.1;
    nbel=m-1;
    ls=dX*((th1-th)/(pi/8)*(1-lsf)+lsf);
    Sna = ls*(1-ra^nbel)/(1-ra)-d12;
       
    while abs(err)>eps
       
       Snb = ls*(1-rb^nbel)/(1-rb)-d12;
       rc = rb-(rb-ra)*Snb/(Snb-Sna);
       Snc = ls*(1-rc^nbel)/(1-rc)-d12;
       err=Snc;
       ra=rb;
       rb=rc;
       Sna=Snb;       
    end
    rbel=rc;
    largestbelow=ls*rbel^nbel;
    below=largestbelow/dX;
    if below>mf
        below
        mf
        pause
    end
   
    
    X0(:,1)=X1;
    X0(:,m)=X2;

    for i=2:nbel
       n=i-2;
       X0(:,i)=X0(:,i-1)+ls*rbel^n*tn;
    end   

    
elseif (th-th2)<=pi/8 && (th-th2)>=0 && lsf<1
   %'above'
    err=1;
    ra=1.05;
    rb=1.1;
    nabv=m-1;
    ls=dX*((th-th2)/(pi/8)*(1-lsf)+lsf);
    Sna = ls*(1-ra^nabv)/(1-ra)-d12;
       
    while abs(err)>eps
       
       Snb = ls*(1-rb^nabv)/(1-rb)-d12;
       rc = rb-(rb-ra)*Snb/(Snb-Sna);
       Snc = ls*(1-rc^nabv)/(1-rc)-d12;
       err=Snc;
       ra=rb;
       rb=rc;
       Sna=Snb;       
    end
    rabv=rc;

    largestafter=ls*rabv^nabv;
    above=largestafter/dX;
    if above>mf
        above
        mf
        pause
    end
    
    X0(:,1)=X1;
    X0(:,m)=X2;

    for i=nabv:-1:2
       n=-i+nabv;
       X0(:,i)=X0(:,i+1)-ls*rabv^n*tn;
    end
          %plot(X0(1,:),X0(2,:),'bx')
    %pause
   
else
    for i=m-1:-1:2
        X0(:,i)=X1+dist(i)*(X2-X1);
    end

end  
    





























 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
m1=m-1; 
x=zeros(2,m); 
x(:,1)=X1; 
dX=(X2-X1)/m1; 
dm=0.*x(1,:);

%
% Apply distribution over X1 to X2
X0=zeros(2,m1+1);
for i=1:m
        X0(:,i)=X1+dist(i)*(X2-X1);
end

    %plot(X0(1,:),X0(2,:),'bx')
    %pause
%}

%if connected<1
for k=1:m1; k1=k+1; 

%   Uniform:    
%   X0=X1+k*dX;
%   [x(:,k1),dm(k1)]=...
%   xy_proj(X0,nn,c1,c2,X,Y,cell,cell_type,cell2p,npc,P,wire,nevtx);

%   Geometric

   [x(:,k1),dm(k1)]=...
   xy_proj(X0(:,k1),nn,c1,c2,X,Y,cell,cell_type,cell2p,npc,P,wire,nevtx);  

end;
%end


dm(1)=dm(2); dm(m)=dm(m1);
xedge=x(1,:); yedge=x(2,:); 


%plot(xedge,yedge,'bx');
%pause
%plot(xedge,yedge,'yx');
%pause

%}

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Projector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [x,dm]=...
    xy_proj(X0,nn,c1,c2,X,Y,cell,cell_type,cell2p,npc,P,wire,nevtx,T) %16



% th    = wire(1);  % Angle of wire
% dw    = wire(2);  % wire diameter
% rf    = wire(3);  % fillet ratio
R     = wire(4);  % Pin radius
% P2    = P/2    ;  % Pitch / 2

npc1 = npc(c1); pc1 = cell2p(c1,1:npc1);
npc2 = npc(c2); pc2 = cell2p(c2,1:npc2);
pcc  = intersect(pc1,pc2);   % Pins that are common to both cells
npcc = length(pcc);


h = .001;

x0=X0; x1=x0+h*nn; 
PP=[x0,x1]';




z1=1; iter=0;
while abs(z1) > 100*eps && iter < 30; iter=iter+1;
   dei=fdedge(PP,X(nevtx+1:nevtx+6),Y(nevtx+1:nevtx+6));
   
   
   if npcc > 1; pcc1=pcc(2:npcc);
     %plot(PP(:,1),PP(:,2),'rx'); %pause
     %plot(X(pcc1),Y(pcc1),'m*');
%     pcc1,pause
      dpi=fdpins(PP,X(pcc1),Y(pcc1),R);
      dw=fdwire(PP,X(pcc1),Y(pcc1),wire,1);
      dp=min(dpi,dw); dei=min(dei,dp); 
   end;

   pcc1 = zeros(1,1); pcc1=pcc(1);
   dpi=fdpins(PP,X(pcc1),Y(pcc1),R);     % subtract lead pin distance
   dw=fdwire(PP,X(pcc1),Y(pcc1),wire,1);
   dp=min(dpi,dw) ;
   dm = min(dp,dei) ;
   de=dp-dei ;
        %plot(X(pcc1),Y(pcc1),'r*');
  %         plot(PP(:,1),PP(:,2),'rx'); %pause
  %   plot(X(pcc(1)),Y(pcc(1)),'m*');
  %   pause
   

   if iter==1; z0=de(1);z1=de(2); else z0=z1; z1=de(1); end;

   dz = z1-z0;
   dx = x1-x0;
   x0 = x1;   
   if z1~=z0;
      x1 = x0 - ( dx./dz )*z1;
      %iter
      %x0,z0
      %x1,z1
      
    else
      iter
      x0,z0
      x1,z1
      'oops'
      pause
   end;
   PP=x1';
end;
x=x0;

theta=atan((x(2)-Y(pcc1))/(x(1)-X(pcc1)));

%plot(x(1)-dpi*cos(theta),x(2)-dpi*sin(theta),'g*')

%plot(x(1)-dw*cos(theta),x(2)-dw*sin(theta),'r*')
%plot(x(1)-dei*cos(theta),x(2)-dei*sin(theta),'kx')
%plot(x(1),x(2),'bx')
%pause


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Return m points marching toward pin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [xedge,yedge]=xy_pin(c1,c4,e,pin,Xcl,Ycl,X,Y,wire,m,T) %17


th    = wire(1);  % Angle of wire
% dw    = wire(2);  % wire diameter
% rf    = wire(3);  % fillet ratio
% R     = wire(4);  % Pin radius

xedge=zeros(m,1);yedge=xedge;xedge(1)=Xcl(c4);yedge(1)=Ycl(c4);

if pin > 0;  % We have an actual pin
    
    
 v1 = [ Xcl(c4)-X(pin) ; Ycl(c4)-Y(pin) ]; v2 = [ cos(th) ; sin(th) ];
 if v1'*v2 > 0;
     
    if T>0
        load wire_quad_mesh_flat; x=xx; y=yy; t=tt; s=ssl; 
    else
        load wire_quad_mesh; x=xx; y=yy; t=tt; s=ssl; 
    end

  [s0,t0,x0,y0]=tmpl_newt(Xcl(c4),Ycl(c4),X(pin),Y(pin),th,xx,yy,tt,ssl);
  m1=m-1; s=s0 + (1-s0)*(0:m1)/m1;
  
  [xedge,yedge]=tmpl_spln(x0,y0,t0,s,X(pin),Y(pin),th,xx,yy,tt,ssl);
  plot(xedge,yedge,'b-'); %pause
 else
  [xedge,yedge]=xy_pin_pr(c1,c4,e,pin,Xcl,Ycl,X,Y,wire,m);
  plot(xedge,yedge,'m-'); %pause
 end;
 
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x]=pin_proj(x0,X,Y,pin,wire) %18

%     Project x0 onto pin+wire

% th    = wire(1);  % Angle of wire
% dw    = wire(2);  % wire diameter
% rf    = wire(3);  % fillet ratio
% R     = wire(4);  % Pin radius


x1=x0;
h=fdpnwr(x1,X(pin),Y(pin),wire); grd=pin_grad(x1,X(pin),Y(pin),wire);
iter=0;
while abs(h) > 100*eps && iter < 30; iter=iter+1;
   x0=x1; x1=x0-h*grd/norm(grd); h=fdpnwr(x1,X(pin),Y(pin),wire);
%  [x1, iter, h]
end;
x=x1;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Take a step of size h in the gradient direction at x0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [x1,d]=pin_grad_step(x0,h,X,Y,pin,wire) %19


grd=pin_grad(x0,X(pin),Y(pin),wire);
x1 = x0+h*grd/norm(grd); d=fdpnwr(x1,X(pin),Y(pin),wire);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Take a step of size h in the pin tangent direction at x0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [x2,x1]=pin_tan_step(x0,h,X,Y,pin,wire) %20


grd=pin_grad(x0,X(pin),Y(pin),wire);
tng=0*grd;tng(1)=-grd(2);tng(2)=grd(1);
x1=x0+h*tng/norm(tng); x2=pin_proj(x1,X,Y,pin,wire);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,R,rw,rf,Db,Xb,thw,thb,fix]=gwire(wire) %21


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
%{
function plt_tmpl2(Xp,Yp,thp) %22

%
%   Load template generated by fem_circ/emesh(); spline; and plot.
%

load wire_quad_mesh; x=xx; y=yy; t=tt; s=ssl; 
[P,R,rw,rf,Db,Xb,thw,thb,fix]=gwire(wire_quad_mesh);

cp=cos(thp);sp=sin(thp);

%
%   Arc-length
%
Ns=length(s); Nt=length(t);

nj=10; nt=Nt;
ns=length(s); at=zeros(ns,nt);
for i=1:ns; 
  for j=2:nt;  % ARCLENGTH
    t0=t(j-1);t1=t(j); dt=(t1-t0)/nj;tt=t0+dt*(0:nj); % nt = LINE RESOLUTION
    xa=spline(t,x(i,:),tt); ya=spline(t,y(i,:),tt); arc=0.;
    for k=2:nj+1; arc=arc+norm([xa(k)-xa(k-1);ya(k)-ya(k-1)]); end;
    at(i,j)=at(i,j-1)+arc;
  end;
end;

%
%  Reparameterize x and y in terms of uniform arclength;
%

%
%   nt lines at constant t;
%

nt=40; ns=length(s); xt=zeros(ns,nt+1);yt=xt; % nt = NUMBER OF LINES
for i=1:ns; 
  t0=min(at(i,:));t1=max(at(i,:)); dt=(t1-t0)/nt; tt = t0 + dt*(0:nt);
  xt(i,:)=spline(at(i,:),x(i,:),tt); 
  yt(i,:)=spline(at(i,:),y(i,:),tt); 
end;

ns=20; s0=min(s);s1=max(s); ds=(s1-s0)/ns; ssl=s0+ds*(0:ns); % ns = LINE RESOLUTION
for j=1:nt+1;
   xs=spline(s,xt(:,j),ssl);
   ys=spline(s,yt(:,j),ssl);
   xx=cp*xs-sp*ys + Xp; yy=sp*xs+cp*ys + Yp; plot(xx,yy,'r-');
end;

%
%   ns lines at constant r;
%
ns=09; s0=min(s);s1=max(s); ds=(s1-s0)/ns; ssl = s0 + ds*(0:ns); % ns = NUMBER OF LINES
nt=length(t); xs=zeros(ns+1,nt); ys=xs;
for j=1:nt; xs(:,j)=spline(s,x(:,j),ssl); ys(:,j)=spline(s,y(:,j),ssl); end;

nt=240; t0=min(t);t1=max(t); dt=(t1-t0)/nt; tt = t0 + dt*(0:nt);% nt = LINE RESOLUTION
for i=1:ns+1;
  xt=spline(t,xs(i,:),tt);
  yt=spline(t,ys(i,:),tt);
  xx=cp*xt-sp*yt + Xp; yy=sp*xt+cp*yt + Yp; plot(xx,yy,'r-');
end;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function plt_tmpl(Xp,Yp,thp) %23


%
%   Load template generated by fem_circ/emesh(); spline; and plot.
%

load wire_quad_mesh;

x=xx; y=yy; t=tt; s=ssl; 

[P,R,rw,rf,Db,Xb,thw,thb,fix]=gwire(wire_quad_mesh);


%
%   Reconstruct high-resolution spline
%
%
%   nt lines at constant t;
%

cp=cos(thp);sp=sin(thp);

nt=40; t0=min(t);t1=max(t); dt=(t1-t0)/nt; tt = t0 + dt*(0:nt); % nt = NUMBER OF LINES
ns=length(s); xt=zeros(ns,nt+1);yt=xt;
for i=1:ns; xt(i,:)=spline(t,x(i,:),tt); yt(i,:)=spline(t,y(i,:),tt); end;

ns=10; s0=min(s);s1=max(s); ds=(s1-s0)/ns; ssl = s0 + ds*(0:ns); % ns = LINE RESOLUTION
for j=1:nt+1;
   xs=spline(s,xt(:,j),ssl);
   ys=spline(s,yt(:,j),ssl);
   xx=cp*xs-sp*ys + Xp; yy=sp*xs+cp*ys + Yp; plot(xx,yy,'r-');
end;

%
%   ns lines at constant r;
%
ns=09; s0=min(s);s1=max(s); ds=(s1-s0)/ns; ssl = s0 + ds*(0:ns); % ns = NUMBER OF LINES
nt=length(t); xs=zeros(ns+1,nt); ys=xs;
for j=1:nt; xs(:,j)=spline(s,x(:,j),ssl); ys(:,j)=spline(s,y(:,j),ssl); end;

nt=240; t0=min(t);t1=max(t); dt=(t1-t0)/nt; tt = t0 + dt*(0:nt);% nt = LINE RESOLUTION
for i=1:ns+1;
  xt=spline(t,xs(i,:),tt);
  yt=spline(t,ys(i,:),tt);
  xx=cp*xt-sp*yt + Xp; yy=sp*xt+cp*yt + Yp; plot(xx,yy,'r-');
end;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [s0,t0,x0,y0]= tmpl_newt(X0,Y0,Xp,Yp,thp,xx,yy,t,s) %24


% Find (s0,t0) s.t. xx(s0,t0) yy(s0,t0) = (X0,Y0);
%% plot(X0,Y0,'cx'); plot(X0,Y0,'co'); plot(X0,Y0,'c*');

%tic
Ns=length(s); Nt=length(t);
cp=cos(thp);sp=sin(thp); xr=cp*xx-sp*yy; yr=sp*xx+cp*yy;
x=xr - (X0-Xp); y=yr - (Y0-Yp);

% Find (s0,t0) s.t. x(s0,t0) y(s0,t0) = (0,0);
d=x.*x+y.*y; d=sqrt(d); d=reshape(d,Ns*Nt,1); [dm,im]=min(d);
i0 = mod(im,Ns); if i0==0; i0=Ns; end; j0=floor((im-1)/Ns)+1;

S1=[s(i0),t(j0)]'; F=[1,1]; iter=0;  % Begin Newton search
%tol = 100*eps;


while norm(F) > 100*eps && iter < 100; iter=iter+1;
   S0=S1;
%  norm(F), iter   %   , toc
   [F,J] = fjac_sp2d(x,y,s,t,S0);
   %** 
   if norm(F) > 100*eps; S1=S0 + J\F; end;
   % If its not converging switch to SOR
   if norm(F) > 100*eps && iter>30; S1=S1-.2*(J\F);end;
   
   x0=F(1)+X0; y0=F(2)+Y0; % plot(x0,y0,'y.',x0,y0,'kx'); pause(.1)
   
end;


s0=S0(1); t0=S0(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [F,J] = fjac_sp2d(xr,yr,s,t,S0) %25


S=S0(1); T=S0(2);

%ns=length(s); nt=length(t);
s0=min(s);s1=max(s); ds=.0001*(s0-s1); t0=min(t);t1=max(t); dt=.0001*(t0-t1);
%s0=min(s);s1=max(s); ds=.00 1*(s0-s1); t0=min(t);t1=max(t); dt=.001*(t0-t1);

%
%  Evaluate ns splines at t0
%
xt0=spline(t,xr,T);yt0=spline(t,yr,T);
ssl=[ S-ds ; S ; S+ds ]; ssl=max(s0,ssl); ssl=min(s1,ssl); ds=max(ssl)-min(ssl);
%%   smin = max(s0,S-ds);
%%   smax = min(s1,S+ds);
%%   smid = (smin+smax)/2.;
%%   ssl=[ smin ; smid; smax ]; ds=max(ssl)-min(ssl);


dxs=spline(s,xt0,ssl); 
dys=spline(s,yt0,ssl);

% 
%  Evaluate nt splines at s0
%
xs0=spline(s,xr',S)'; ys0=spline(s,yr',S)';
tt=[ T-dt ; T ; T+dt ]; tt=max(t0,tt); tt=min(t1,tt); dt=max(tt)-min(tt);
dxt=spline(t,xs0,tt); dyt=spline(t,ys0,tt);

F = .5*[ dxs(2)+dxt(2); dys(2)+dyt(2) ];


ds = max(ssl)-min(ssl); dt = max(tt)-min(tt);
dFxds = (dxs(3)-dxs(1))/ds; dFyds = (dys(3)-dys(1))/ds;
dFxdt = (dxt(3)-dxt(1))/dt; dFydt = (dyt(3)-dyt(1))/dt;

J = [ dFxds , dFxdt ; dFyds , dFydt ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xedge,yedge]=xy_pin_pr(c1,c4,e,pin,Xcl,Ycl,X,Y,wire,m) %26


% th    = wire(1);  % Angle of wire
% dw    = wire(2);  % wire diameter
% rf    = wire(3);  % fillet ratio
% R     = wire(4);  % Pin radius

X0=[Xcl(c4),Ycl(c4)];
xedge=zeros(m,1);yedge=xedge;xedge(1)=X0(1);yedge(1)=X0(2);

if pin > 0;  % We have an actual pin
   X1=pin_proj(X0,X,Y,pin,wire); dx=(X1-X0)/(m-1); h=norm(dx);
   err=1; iter=0;
   while err > 100*eps && iter < 30; iter=iter+1;
      x1=X0;
      for k=2:m
        x0=x1; [x1,d]=pin_grad_step(x0,-h,X,Y,pin,wire);
        xedge(k)=x1(1); yedge(k)=x1(2);
      end;
      plot(xedge,yedge,'m-');  %HERE
      
      err=abs(d); h = h + d/m;
%     [x1,iter,h,d]
%     plot(x1(1),x1(2),'bd'), pause(.1); pause
   end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xedge,yedge]=tmpl_spln(x0,y0,t0,s,X,Y,th,xx,yy,tt,ssl) %27


Ns=length(ssl); Nt=length(tt); cp=cos(th);sp=sin(th);


%  Evaluate ns splines at t0

xt0=zeros(Ns,1); yt0=xt0;


xt0(:)=spline(tt,xx(:,:),t0); 
yt0(:)=spline(tt,yy(:,:),t0);

%clf
%axis equal
%hold on
%plot(xx,yy,'b-')


%xt0(1:Ns)=spline(tt,xx(1:Ns,:),t0); yt0(1:Ns)=spline(tt,yy(1:Ns,:),t0); 


xedge=spline(ssl,xt0,s); yedge = spline(ssl,yt0,s);

%plot(xedge,yedge,'mx')
%pause

xr=cp*xedge-sp*yedge; yr=sp*xedge+cp*yedge;
xedge=xr+X; yedge=yr+Y;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xedge,yedge,tp]=xy_pin2(pin,X0,Y0,X,Y,wire,m,dist,xx,yy,tt,ssl) %28


th    = wire(1);  % Angle of wire
% dw    = wire(2);  % wire diameter
% rf    = wire(3);  % fillet ratio
% R     = wire(4);  % Pin radius

xedge=zeros(m,1);yedge=xedge;xedge(1)=X0;yedge(1)=Y0;

if pin > 0;  % We have an actual pin

  v1 = [ X0-X(pin) ; Y0-Y(pin) ]; v2 = [ cos(th) ; sin(th) ];
  if v1'*v2 > 0;
    %load wire_quad_mesh; x=xx; y=yy; t=tt; s=ssl; 

    [s0,t0,x0,y0]= tmpl_newt(X0,Y0,X(pin),Y(pin),th,xx,yy,tt,ssl);


    
    %m1=m-1; s=s0 + (1-s0)*(0:m1)/m1;
       % different distribution
    s=s0+(1-s0)*dist;
    
    [xedge,yedge]=tmpl_spln(x0,y0,t0,s,X(pin),Y(pin),th,xx,yy,tt,ssl);
    

    
    %[xedge,yedge,sa]=spline_block_arc1(xedge,yedge); % FIXED ARC SPACING
    tp=1;
  else
    [xedge,yedge]=xy_pin_pr2(pin,X0,Y0,X,Y,wire,m,dist);
   % plot(X0,Y0,'r*')
    tp=2;
   % plot(xedge,yedge,'bx')
   % pause
  end;

  
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xedge,yedge]=xy_pin_pr2(pin,x0,y0,X,Y,wire,m,dist) %29


% th    = wire(1);  % Angle of wire
% dw    = wire(2);  % wire diameter
% rf    = wire(3);  % fillet ratio
% R     = wire(4);  % Pin radius

X0=[x0,y0];
xedge=zeros(m,1);yedge=xedge;xedge(1)=X0(1);yedge(1)=X0(2);

if pin > 0;  % We have an actual pin
   X1=pin_proj(X0,X,Y,pin,wire); dx=(X1-X0)/(m-1); %h=norm(dx);
   err=1; iter=0;
   
   Xd=zeros(2,m);
   %hd=zeros(1,m);
   
   for i=1:m
        Xd(:,i)=X0+dist(i)*(X1-X0);
   end
   xedge=Xd(1,:)';
   yedge=Xd(2,:)';

%{
   while err > 100*eps && iter < 100; iter=iter+1;
      
      x1=X0;
      for k=2:m
        h=norm(Xd(:,k)-Xd(:,k-1));
        x0=x1; [x1,d]=pin_grad_step(x0,-h,X,Y,pin,wire);
        xedge(k)=x1(1); yedge(k)=x1(2);
      end;
     %plot(xedge,yedge,'mx');
      err=abs(d); h = h + d/m;
%     [err,iter]
%     [x1,iter,h,d]
%     plot(x1(1),x1(2),'bd'), pause(.1);% pause
   
    end;
%}

%else;
%   'no pin'
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spline_block(X,Y) %30


    m=size(X,1); n=size(X,2); s=0.1;

    xx=spline(1:n,X,1:s:n); yy=spline(1:n,Y,1:s:n);
%   for i=1:m; plot(xx(i,:),yy(i,:),'m-'); end;
    for i=1:m-1:m; plot(xx(i,:),yy(i,:),'m-'); end;

    xx=spline(1:m,X',1:s:m); yy=spline(1:m,Y',1:s:m);
%   for j=1:n; plot(xx(j,:),yy(j,:),'m-'); end; pause (.01);
    for j=1:n-1:n; plot(xx(j,:),yy(j,:),'m-'); end; %pause (.01);
 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xa,Ya,sa]=spline_block_arc1(X,Y) %31


n=length(X); ds=0.1;

s=0:(n-1); ssl=0:ds:n-1; ns=length(ssl); aa=arclength(s,X,Y,ssl); 
aa=(n-1)*aa./aa(ns); 
a=s; sa = spline(aa,ssl,a); % Uniform sample in archlength;
Xa=spline(s,X,sa); Ya=spline(s,Y,sa); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xa,Ya,sa]=spline_block_arc(X,Y) %32


m=size(X,1); n=size(X,2); ds=0.1;

s=0:(n-1); ssl=0:ds:n-1; ns=length(ssl); aa=arclength(s,X(m,:),Y(m,:),ssl); 
aa=(n-1)*aa./aa(ns); 
a=s; sa = spline(aa,ssl,a); % Uniform sample in archlength;
Xa=spline(s,X,sa); Ya=spline(s,Y,sa); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aa=arclength(s,x,y,ssl) %33


xx=spline(s,x,ssl); yy=spline(s,y,ssl);
dx=diff(xx);dy=diff(yy);aa=dx.*dx+dy.*dy;aa=sqrt(aa);
[n1,n2]=size(aa);
if n1>n2; aa=[0;aa]; else aa=[0,aa]; end; aa=cumsum(aa);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pst]=pt_on_line(p,x1,x2,count) %34

th  = (x2-x1)/norm(x2-x1); pst =  x1 + ((p-x1)'*th)*th;

%% if count==1;
%%    plot([pst(1),p(1)],[pst(2),p(2)],'r-');
%%    plot(pst(1),pst(2),'ro');
%%    plot(p(1),p(2),'r*');
%% end; if count==2;
%%    plot([pst(1),p(1)],[pst(2),p(2)],'y-');
%%    plot(pst(1),pst(2),'yo');
%%    plot(p(1),p(2),'y*');
%% end; if count==3;
%%    plot([pst(1),p(1)],[pst(2),p(2)],'k-');
%%    plot(pst(1),pst(2),'ko');
%%    plot(p(1),p(2),'k*');
%% end; if count==4;
%%    plot([pst(1),p(1)],[pst(2),p(2)],'b-');
%%    plot(pst(1),pst(2),'bo');
%%    plot(p(1),p(2),'b*');
%% end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ...
[nel,nelp,nele,npts,ve,ce,be,pe,qe,te,Xcl,Ycl,neigh_e,neigh_s]=...
  gen_block(ncell,npin,ne,ncp,p2cell,X,Y,Xcl,Ycl,R,ntcl,necl,ecell,ccell,T) %35


% Input:
%
%    X,Y     --- cell vertex coordinates (includes pin centers for triangular cells)
%    Xcl,Ycl --- cell center coordinates
%
% Output:
%
%    Xcl,Ycl --- cell center coordinates + edge & corner nodes that make basic
%                SEM building blocks
%
% GENERATE CELL-BASED MESH.

nel=0; 
ve=zeros(4,12*npin);    % vertex pointers for each element, e=1:nel
ce=ve;                  % curve side info for each element
be=ve;                  % boundary conditions for each element
pe=zeros(1,12*npin);    % pin id for each element
qe=zeros(1,12*npin);    % pin id for each element


npts=ncell;
for k=1:npin; ncl=ncp(k); Xp=X(k); Yp=Y(k); xp=[Xp,Yp]';

  % BUILD PIN-BASED ELEMENTS

  % sort cell nodes by angle
  x=Xcl(p2cell(k,1:ncl)); y=Ycl(p2cell(k,1:ncl)); t=atan2(y-Yp,x-Xp);
  
  
  for i=1:ncl; 
      if t(i)<0; t(i)=t(i)+2*pi;end; 
      
      
      
  end;
  
  
  [s,I]=sort(t); x=x(I); y=y(I); p2cell(k,1:ncl)=p2cell(k,I);% plot(x,y,'k-');

  if k==12
      t
  end
  

  i=ncl;
  x4=[Xcl(p2cell(k,i)),Ycl(p2cell(k,i))]';
  r4=(x4-xp)/norm(x4-xp); x1=xp+R*r4;%EB notation
  npts=npts+1; Xcl(npts)=x1(1); Ycl(npts)=x1(2); npts1=npts; 
  % check for speed
  for i=1:ncl; i1=i-1; if i1==0; i1=ncl; end; nel=nel+1;
    x2=x1; x3=x4;
    x4=[Xcl(p2cell(k,i)),Ycl(p2cell(k,i))]';
    r4=(x4-xp)/norm(x4-xp); x1=xp+R*r4;%non-sym. notation
    if i < ncl;
      npts=npts+1; Xcl(npts)=x1(1); Ycl(npts)=x1(2);
      ve(:,nel)= [npts,npts-1,p2cell(k,i1),p2cell(k,i)]';
    else
      ve(:,nel)= [npts1,npts,p2cell(k,i1),p2cell(k,i)]';
    end;
    ce(1,nel)= -R; pe(nel)=k; qe(nel)=k; te(nel)=1;  % Curve, Pin, Type
    txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'r-');hold on;
  end;
end;   nelp = nel; nptsp = npts;


%    GENERATE EDGE ELEMENTS

Re=2.5*R; % Just a guess for edge radius
% Check for speed
for iedge=1:6; 
    for ee=1:ne-1;
       if ee==1; icl=1+(ne-1)*(iedge-1);
         xc1=[X(ecell(1,icl)),Y(ecell(1,icl))]'; 
         xc2=[X(ecell(2,icl)),Y(ecell(2,icl))]';
         xe2=pt_on_line([Xcl(ntcl+icl),Ycl(ntcl+icl)]',xc1,xc2,3);
         npts=npts+1; Xcl(npts)=xe2(1); Ycl(npts)=xe2(2);
       end;
       if ee < ne-1;
         icl=icl+1;
         xc1=[X(ecell(1,icl)),Y(ecell(1,icl))]'; 
         xc2=[X(ecell(2,icl)),Y(ecell(2,icl))]';
         xe2=pt_on_line([Xcl(ntcl+icl),Ycl(ntcl+icl)]',xc1,xc2,4);
         npts=npts+1; Xcl(npts)=xe2(1); Ycl(npts)=xe2(2);

         nel=nel+1; 
         ve(:,nel)=[npts-1,npts,ntcl+icl,ntcl+icl-1]'; ce(3,nel)= -Re; 
         
         
         pe(nel)=-((ne-2)*(iedge-1)+ee); qe(nel)=-((ne-2)*(iedge-1)+ee)-6; te(nel)=2;
         txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'r-');
       end;
    end;
end;  nele = nel-nelp; nptse=npts-nptsp;


%    GENERATE CORNER ELEMENTS

%
%   For corner and edge blocks, pe < 0, and |pe| points to the 
%   subassembly edge associated with side 1 of block e.
%

Rc=2.5*R; % Just a guess for corner radius

if ne==2; 
  for icrn=1:6; 
    npts=npts+1; Xcl(npts)=X(ccell(2,icrn)); Ycl(npts)=Y(ccell(2,icrn));
    if icrn==1;
      nel=nel+1;
      ve(:,nel)=[nptsp+6,npts,ntcl+necl+icrn,ntcl+necl]';ce(3,nel)= -Rc;
      pe(nel)=-6; qe(nel)=-1; te(nel)=3;
      txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'y-');
      nel=nel+1; 
      ve(:,nel)= [npts,nptsp+1,ntcl+1,ntcl+necl+icrn]'; ce(3,nel)= -Rc;
      pe(nel)=-1; qe(nel)=-1; te(nel)=3;
      txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'y-');
    else
      nel=nel+1; 
      ve(:,nel)=[nptsp+icrn-1,npts,ntcl+necl+icrn,ntcl+icrn-1]'; ce(3,nel)= -Rc;
      pe(nel)=-(icrn-1); qe(nel)=-icrn; te(nel)=3;
      txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'y-');
      nel=nel+1; 
      ve(:,nel)=[npts,nptsp+icrn,ntcl+icrn,ntcl+necl+icrn]';ce(3,nel)=-Rc;
      pe(nel)=-icrn; qe(nel)=-icrn; te(nel)=3;
      txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'y-');
    end;
  end;
else 
  for icrn=1:6; 
    npts=npts+1; Xcl(npts)=X(ccell(2,icrn)); Ycl(npts)=Y(ccell(2,icrn));

    if icrn==1;
     nel=nel+1;
     ve(:,nel)=[ve(2,nelp+nele),npts,ntcl+necl+icrn,ve(3,nelp+nele)]';
     ce(3,nel)= -Rc; pe(nel)=-6; qe(nel)=-1; te(nel)=3;
     txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'y-');
     nel=nel+1; 
     ve(:,nel)= [npts,ve(1,nelp+1),ve(4,nelp+1),ntcl+necl+icrn]';
     ce(3,nel)= -Rc; pe(nel)=-1; qe(nel)=-1; te(nel)=3;
     txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'y-');
    else
     el_edge = nelp+(ne-2)*(icrn-1);
     nel=nel+1; 
     ve(:,nel)= [ve(2,el_edge),npts,ntcl+necl+icrn,ve(3,el_edge)]';
     ce(3,nel)= -Rc; pe(nel)=-icrn; qe(nel)=-icrn; te(nel)=3;
     txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'y-');
     nel=nel+1; 
     ve(:,nel)= [npts,ve(1,el_edge+1),ve(4,el_edge+1),ntcl+necl+icrn]';
     icrn1=icrn+1; if icrn1>6; icrn1=1; end;
     ce(3,nel)= -Rc; pe(nel)=-icrn1; qe(nel)=-icrn; te(nel)=3;
     txt_elem(Xcl,Ycl,ve(:,nel),ce(:,nel),sprintf('%d',pe(nel)),'y-');
    end;
  end;
end; 
%nelc = 12;


[neigh_s,neigh_e]=get_el_neighbors(ve,nel);

for ie=1:nel; 
    for is=1:4;      %  Clean up curve sides
        je = neigh_e(is,ie); js = neigh_s(is,ie);
        if je>0; if abs(ce(is,ie))>0; ce(js,je)=-ce(is,ie); end; end;
    end; 
end;

for ie=1:nel; %  Clean up domain Types (for proper element decomposition)
    je = neigh_e(3,ie); 
    te(ie) = max(te(ie),te(je)); te(je)=te(ie);
    txt_elem(Xcl,Ycl,ve(:,ie),ce(:,ie),sprintf('%d',te(ie)),'y-');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nvtx,nevtx,ntcl,tcell,necl,ecell,nccl,ccell,X,Y]=...
    get_cells(ne,npin,km,xm,ym,X,Y,D,P) % Construct cells connecting pins %36



PD = P/D; R=D/2; p2flat = R + D*(PD-1);
nr = 2*ne-1;                     % # pin rows and max in any row
nc = ne;

tcell = zeros(3,2*nr*nr);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Construct triangular cells connecting pins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
k=0; ntcl=0;  % ntcl = # ncells
for ir=1:nr;
  for jc=1:nc;
    k=k+1;
    if jc < nc && ir < ne;
      ntcl=ntcl+1; tcell(:,ntcl)=[k,km(ir+1,jc+1),km(ir+1,jc  )]';
      ntcl=ntcl+1; tcell(:,ntcl)=[k,k+1          ,km(ir+1,jc+1)]';
    elseif jc == nc && ir < ne;
      ntcl=ntcl+1; tcell(:,ntcl)=[k,km(ir+1,jc+1),km(ir+1,jc  )]';
    elseif jc == 1 && ir < nr;
      ntcl=ntcl+1; tcell(:,ntcl)=[k,k+1          ,km(ir+1,jc  )]';
    elseif jc < nc && ir < nr;
      ntcl=ntcl+1; tcell(:,ntcl)=[k,km(ir+1,jc  ),km(ir+1,jc-1)]';
      ntcl=ntcl+1; tcell(:,ntcl)=[k,k+1          ,km(ir+1,jc  )]';
    end;
  end;
  if ir<ne; nc=nc+1; elseif ir<nr; nc=nc-1; end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Define edge cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th=pi*(0:5)/3; xn=sin(th);yn=-cos(th);

nvtx=npin; 
necl=0; ecell=zeros(4,6*(ne-1));

ke=1; ir=1;              % Bottom row
for jc=1:nc;
  nvtx=nvtx+1; X(nvtx)=xm(ir,jc)+p2flat*xn(ke); Y(nvtx)=ym(ir,jc)+p2flat*yn(ke);
  if jc>1; necl=necl+1;ecell(:,necl)=[nvtx-1,nvtx,km(ir,jc),km(ir,jc-1)]';end;
end;

ke=2; 
for ir=1:ne;       % Lower right edge
  nvtx=nvtx+1; X(nvtx)=xm(ir,nc)+p2flat*xn(ke); Y(nvtx)=ym(ir,nc)+p2flat*yn(ke);
  if ir>1; necl=necl+1;ecell(:,necl)=[nvtx-1,nvtx,km(ir,nc),km(ir-1,nc-1)]';end;
  if ir<ne; nc=nc+1; end;
end;

ke=3; 
for ir=ne:nr;      % Upper right edge
  nvtx=nvtx+1; X(nvtx)=xm(ir,nc)+p2flat*xn(ke); Y(nvtx)=ym(ir,nc)+p2flat*yn(ke);
  if ir>ne;necl=necl+1;ecell(:,necl)=[nvtx-1,nvtx,km(ir,nc),km(ir-1,nc+1)]';end;
  if ir<nr; nc=nc-1; end;
end;

ke=4; ir=nr;             % Top row
for jc=nc:-1:1;
  nvtx=nvtx+1; X(nvtx)=xm(ir,jc)+p2flat*xn(ke); Y(nvtx)=ym(ir,jc)+p2flat*yn(ke);
  if jc<nc;necl=necl+1;ecell(:,necl)=[nvtx-1,nvtx,km(ir,jc),km(ir,jc+1)]';end;
end;

ke=5; 
for ir=nr:-1:ne;   % Upper left edge
  nvtx=nvtx+1; X(nvtx)=xm(ir,1)+p2flat*xn(ke); Y(nvtx)=ym(ir,1)+p2flat*yn(ke);
  if ir<nr;necl=necl+1;ecell(:,necl)=[nvtx-1,nvtx,km(ir,1),km(ir+1,1)]';end;
end;

ke=6; 
for ir=ne:-1:1;   % Lower left edge
  nvtx=nvtx+1; X(nvtx)=xm(ir,1)+p2flat*xn(ke); Y(nvtx)=ym(ir,1)+p2flat*yn(ke);
  if ir<ne;necl=necl+1;ecell(:,necl)=[nvtx-1,nvtx,km(ir,1),km(ir+1,1)]';end;
end;
nevtx = nvtx;



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Define corner cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
nccl=6; ccell=zeros(4,nccl);
Rsmbly = (max(Y)-min(Y))/sqrt(3.); th=pi*(-2:3)/3; ct=cos(th);st=sin(th);
for k=1:6;
  c=1+(ne-1)*(k-1); c1=c-1; if k==1; c1=necl; end;
  nvtx=nvtx+1; X(nvtx)=Rsmbly*ct(k); Y(nvtx)=Rsmbly*st(k);
  ccell(:,k)=[ecell(2,c1),nvtx,ecell(1,c),ecell(4,c)]';
end;

%  for k=1:ntcl; txt_cell(X,Y,tcell(:,k),sprintf('%d',k)); hold on; end;
%  for k=1:necl; txt_cell(X,Y,ecell(:,k),sprintf('%d',k+ntcl)); hold on; end;
%  for k=1:nccl; txt_cell(X,Y,ccell(:,k),sprintf('%d',k+ntcl+necl));hold on;end;
%  'CELLS, PAUSE'
%  pause


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Define cell centers ( := maximum pt of signed distance function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xcl,Ycl]=cell_center(npin,X,Y,ntcl,tcell,necl,ecell,nccl,ccell,R) %37


%n=length(X);
%plot(X,Y,'ro')
%'cell_center'
%pause

ncell = ntcl+necl+nccl;
Xcl=zeros(ncell+6*npin,1); Ycl=Xcl;

%                              % INTERIOR CELLS
for k=1:ntcl; Xcl(k)=sum(X(tcell(:,k)))/3; Ycl(k)=sum(Y(tcell(:,k)))/3; end;
% Check for speed
for j=1:necl; k=j+ntcl;        % EDGE CELLS
    x1 = [ X(ecell(1,j)) , Y(ecell(1,j)) ];
    x2 = [ X(ecell(2,j)) , Y(ecell(2,j)) ];
    xp = [ X(ecell(3,j)) , Y(ecell(3,j)) ];
    dt = x2-x1; lt=norm(dt);
    dn = xp-x2; ln=norm(dn);
    an = (ln.^2 + (lt/2).^2-R.^2)/(2*(ln+R));
    Xcl(k)=sum(X(ecell(1:2,j)))/2; Ycl(k)=sum(Y(ecell(1:2,j)))/2;
    Xcl(k)=Xcl(k) + an*dn(1)/ln; Ycl(k)=Ycl(k) + an*dn(2)/ln; 
end;

for j=1:nccl; k=j+ntcl+necl;   % CORNER CELLS
    x2 = [ X(ccell(2,j)) , Y(ccell(2,j)) ];
    x4 = [ X(ccell(4,j)) , Y(ccell(4,j)) ];
    dr = x4-x2; lr=norm(dr);
    l = (lr-R)/(1.+2./sqrt(3.)); s = 2.*l/sqrt(3.);
    Xcl(k)=x2(1)+s*dr(1)/lr; Ycl(k)=x2(2)+s*dr(2)/lr; 
end;

%%plot(Xcl(1:ncell),Ycl(1:ncell),'rx'); pause
