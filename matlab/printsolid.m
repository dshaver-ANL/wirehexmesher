%printsolid
load savesolid9


fopen( 'pxmesh_8-14.out', 'wt' );
fopen( 'pymesh_8-14.out', 'wt' );
fopen( 'pzmesh_8-14.out', 'wt' );
fopen( 'pdata_v5f.out', 'wt' );

tic
dlmwrite('pdata_v5f.out',r_data,'precision', '%12.8f','-append');
toc
tic
%dlmwrite('pxmesh_simple.out',pmeshx,'precision', '%10.6f','-append');
%dlmwrite('pymesh_simple.out',pmeshy,'precision', '%10.6f','-append');
%dlmwrite('pzmesh_simple.out',pmeshz,'precision', '%10.6f','-append');
dlmwrite('pxmesh_8-14.out',pmeshx,'precision', '%12.8f','-append');
toc
dlmwrite('pymesh_8-14.out',pmeshy,'precision', '%12.8f','-append');
toc
dlmwrite('pzmesh_8-14.out',pmeshz,'precision', '%12.8f','-append');
toc
