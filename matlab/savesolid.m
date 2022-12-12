


load savesolid 

%load -mat 61pin_12x3.mat  % pin_blcko, pin_block, type_block, xyzelems, ne, Dw, Df

P=18.875;
D=15.875;



eps_d2=5*1e-6;
PD=P/D;
%ne=3;
Do=1;
Di=5.65/6.55;
Dwo=Dw;%1.7/6.55;
Dwi=0.12;%Dwo-(Do-Di);
%Df=.2/6.55;
%Nx=4;

[npin,X,Y,xm,ym,km]=pin_coors(ne,PD);



npins=max(pin_blcko);
Nx=size(xelems,1)-1;
Col=size(xelems,5);
Lay=size(xelems,7);
nblocks=size(xelems,6)-6*ne; 



tic
for level=1:(Lay-1)
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



level=Lay;
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
                  pmeshz(elemnum,gr,gc,3)=zelems(1,gc,1,1,j,i+block1-1,1);    
                  end;
                end
               
            end
         end

end;
save savesolid2 pmeshx pmeshy pmeshz
