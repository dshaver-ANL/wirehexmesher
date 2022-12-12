  function [to,po,ind] = reorder_p(t,p,list); % points in list ordered last

  nb = max(max(t));
  nl = length(list);

  ip = 1:nb; ip(list) = nb+1:nb+nl;

  [is,ind  ] = sort(ip);   % ind  lists _former_ positions
  [is,rankl] = sort(ind);  % rank gives orignal  ranking
% min(rankl)
% pause

  po = p(ind,:); 
  to = rankl(t);
  

