% one_pin.m

% Takes the information from xyzelems to create a binary file with only the
% information of the first pin

load ess_xyblocks % pin_blcko, pin_block, type_block, xyzelems
% xyzelems(3, GLL R, GLL C, GLL L, Row, Column, elem, Layer)


% The elements that describe the first pin are:
ec=(1);
nec=size(ec,2);

pin=xyzelems(:,:,:,:,:,:,ec,1);

fid = fopen('one_element_l.bin','w');
count = fwrite(fid, pin, 'single', 0, 'l');
count




% Display it
%%{
pin_old=1;
for L=1%:Lay
   for gl=1:Nx+1   
      clf
      hold on
      axis equal
      for e=1:nec
         for C=1:Col
            for R=1:Row
               X=squeeze(pin(1,:,:,gl,R,C,e,L));Y=squeeze(pin(2,:,:,gl,R,C,e,L));
               m=size(X,1);n=size(X,2); s=0.1;
               xx=spline(1:n,X,1:s:n); yy=spline(1:n,Y,1:s:n);
               for i=1:m-1:m; plot(xx(i,:),yy(i,:),'m-'); end;
               xx=spline(1:m,X',1:s:m); yy=spline(1:m,Y',1:s:m);
               for j=1:n-1:n; plot(xx(j,:),yy(j,:),'m-'); end; 

               plot(X,Y,'bx')
            end
         end
      end
      pause
   end
end
%}











