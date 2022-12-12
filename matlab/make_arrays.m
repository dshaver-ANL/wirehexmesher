function [cpin, epin]=make_arrays(ne)

%ne=3; % number of edge pins
npt=3*ne^2-3*ne+1; % total numper of pins
npi=6*ne-5; % number of important pins


% cpin is an array that renumbers the important pins and distinguishes them
% from unimportant [ 1 2 3 4 0 0 5 6 0 7 0 8 9 0 0 10 11 12 13];
cpin=zeros(1,npt);
p=0; % place
v=0; % value
for j=1:ne; p=p+1;v=v+1; cpin(p)=v;end
for k=2:ne-1
   for j=1:k+ne-1
      if j==1||j==k+ne-1;p=p+1;v=v+1;cpin(p)=v;
      else p=p+1;cpin(p)=0;end
   end
end
p=p+1;v=v+1;cpin(p)=v;
for j=1:ne-2; p=p+1;cpin(p)=0; end
p=p+1;v=v+1;cpin(p)=v;
for j=1:ne-2; p=p+1;cpin(p)=0; end
p=p+1;v=v+1;cpin(p)=v;
for k=ne-1:-1:2
   for j=1:k+ne-1
      if j==1||j==k+ne-1; p=p+1;v=v+1;cpin(p)=v;
      else p=p+1;cpin(p)=0; end
   end
end
for j=1:ne; p=p+1;v=v+1;cpin(p)=v; end

% ccw is an array that renumbers the important pins ccw
% [1 2 3 5 8 10 13 12 11 9 6 4 7]
ccw=zeros(1,npi);
ccw(1)=1;                                                       
for i=2:ne; ccw(i)=ccw(i-1)+1; end
for i=ne+1:3*ne-3; ccw(i)=ccw(i-1)+2; end
for i=2*ne-1:3*ne-3; ccw(i)=ccw(i)+1; end
ccw(3*ne-2)=ccw(3*ne-3)+ne;
for i=3*ne-1:4*ne-3; ccw(i)=ccw(i-1)-1; end
for i=4*ne-2:6*ne-6; ccw(i)=ccw(i-1)-2; end
for i=5*ne-4:6*ne-6; ccw(i)=ccw(i)-1; end
ccw(6*ne-5)=3*ne-2;


% epin gives an order for the ouside elements
% [ 1 3 8 13 11 6 2 5  10 12 9 4];
epin=zeros(1,npi-1);
for i=1:6
   epin(i)=ccw((ne-1)*i-(ne-2));
   
   for j=1:ne-2
      epin(8-ne+i*(ne-2)+j)=ccw((ne-1)*i+j-ne+2);   
   end
end


















