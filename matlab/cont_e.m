  function [x,a] = cont_e(t,X,e,u,uv); % single element contour line


  s = u(t(e,:)) ;
  p = X(t(e,:),:);

  k=0;
  a=0;
  for j=1:3;
     j1 = j+1;  if j1>3 , j1=1; end;
     if (s(j)-uv)*(s(j1)-uv) <=0 & s(j) ~= s(j1)
        k = k+1;
        a(k) = (uv-s(j))/(s(j1)-s(j));
        x(k,:) = p(j,:) + a(k)*(p(j1,:)-p(j,:));
     end;
  end;