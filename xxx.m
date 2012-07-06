
numRepeats = 4;

IPTGFiles = [ '1'; '2'; '3' ];
IPTGConcs = { '10µM'; '20µM'; '40µM' };

arabinoseFiles = [ '4'; '5'; '6' ];
arabinoseConcs = { '0.005%'; '0.01%'; '0.02%' };

cheYFile = 'CheY';
cheZFile = 'CheZ';

docroot = '/Nikon/2012/06/25/';

i = 1;
cellCount = zeros( 1, 5 );
for j = 1:numRepeats
    cell_loc =...
        strcat( docroot, IPTGFiles(i), '-', cheYFile, num2str(j), '.tif' );
    cluster_loc =...
        strcat( docroot, IPTGFiles(i), '-', cheZFile, num2str(j), '.tif' );

    ci = ColiImage( cell_loc, cluster_loc, 160 );

    cellCount = cellCount + ci.clusterHistogram();
end

bar( 0:4, cellCount, 'hist' );