function [Ab,Q]=afem(p,t) % Build 2D FEM stiffness matrix

  nt = size(t,1); nl=3*nt;


  y23 = p(t(:,2),2) - p(t(:,3),2);
  y31 = p(t(:,3),2) - p(t(:,1),2);
  y12 = p(t(:,1),2) - p(t(:,2),2);
  x32 = p(t(:,3),1) - p(t(:,2),1);
  x13 = p(t(:,1),1) - p(t(:,3),1);
  x21 = p(t(:,2),1) - p(t(:,1),1);
  area2 = x21.*y31 - y12.*x13 ;
  arean = min(area2);
  aream = max(area2);


  for e=1:nt;
    if area2(e) < 0;
       temp   = t(e,2);
       t(e,2) = t(e,3);
       t(e,3) = temp;
    end;
  end;

  y23 = p(t(:,2),2) - p(t(:,3),2);
  y31 = p(t(:,3),2) - p(t(:,1),2);
  y12 = p(t(:,1),2) - p(t(:,2),2);
  x32 = p(t(:,3),1) - p(t(:,2),1);
  x13 = p(t(:,1),1) - p(t(:,3),1);
  x21 = p(t(:,2),1) - p(t(:,1),1);
  area4 = (0.50./(x21.*y31 - y12.*x13)) ;



  i0 = (0:nt-1)';
  i1 = (1:nt)';
  AL = spalloc(nl,nl,3*nl);
  A1 = zeros(3,3,nt);

  A1(1,1,:) = area4.*( y23.*y23+x32.*x32 );
  A1(1,2,:) = area4.*( y23.*y31+x32.*x13 );
  A1(1,3,:) = area4.*( y23.*y12+x32.*x21 );
  A1(2,1,:) = area4.*( y31.*y23+x13.*x32 );
  A1(2,2,:) = area4.*( y31.*y31+x13.*x13 );
  A1(2,3,:) = area4.*( y31.*y12+x13.*x21 );
  A1(3,1,:) = area4.*( y12.*y23+x21.*x32 );
  A1(3,2,:) = area4.*( y12.*y31+x21.*x13 );
  A1(3,3,:) = area4.*( y12.*y12+x21.*x21 );

  for e=0:nt-1;
    AL(3*e+(1:3),3*e+(1:3)) = A1(:,:,e+1);
  end;

  Q  = sparse(1:nl,reshape(t',nl,1),1);
  Ab = Q'*AL*Q;