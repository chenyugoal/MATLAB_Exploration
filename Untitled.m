Spiketimes = load('SpikeTrain2.txt')

%box = zeros(12,659); % break original data into 12 spike trains artificially
%for i=1:659
%    box(:,i) = Spiketimes((i-1)*12+1:i*12)
%end

histogram(Spiketimes,'Normalization','count','BinWidth',0.008,'DisplayStyle','stairs')