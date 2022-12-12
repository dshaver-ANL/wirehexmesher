function [to]=chk_jac(t,p) % Reorder local triangles if area(e) < 0

  nt = size(t,1);


  y23 = p(t(:,2),2) - p(t(:,3),2);
  y31 = p(t(:,3),2) - p(t(:,1),2);
  y12 = p(t(:,1),2) - p(t(:,2),2);
  x32 = p(t(:,3),1) - p(t(:,2),1);
  x13 = p(t(:,1),1) - p(t(:,3),1);
  x21 = p(t(:,2),1) - p(t(:,1),1);
  area2 = x21.*y31 - y12.*x13 ;
  arean = min(area2);
  aream = max(area2);


  to = t;
  for e=1:nt;
    if area2(e) < 0;
       to(e,2) = t(e,3);
       to(e,3) = t(e,2);
    end;
  end;