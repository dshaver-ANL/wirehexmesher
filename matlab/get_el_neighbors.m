function [neigh_s,neigh_e]=get_el_neighbors(cell,nel)
neigh_s = 0.*cell; neigh_e = neigh_s;

edge = zeros(4*nel,1); neighs_e = edge; neighs_s = edge;

emx = max(max(cell));
k=0; for e=1:nel; for i=1:4; i1=i+1; if i1>4; i1=1; end;
   e1=min(cell(i,e),cell(i1,e));
   e2=max(cell(i,e),cell(i1,e));
   k=k+1;
   edge   (k)=e1 + emx*(e2-1);
   neighs_e(k)=e;
   neighs_s(k)=i;
end;end;
[edge,I] = sort(edge); neighs_e=neighs_e(I); neighs_s=neighs_s(I);

k=1; while k < 4*nel; 
  k1=k+1;
  if edge(k)==edge(k1);
     ie=neighs_e(k ); is=neighs_s(k ); 
     je=neighs_e(k1); js=neighs_s(k1); 
     neigh_e(is,ie) = je; neigh_s(is,ie) = js;
     neigh_e(js,je) = ie; neigh_s(js,je) = is;
     k=k+1;
  end;
  k=k+1;
    
end;

