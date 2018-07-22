function [aligned] = pp_DataSets(ca)

nca = length(ca);
colors=['k','r','b','g','m','c'];

for j=1:nca
  [n,m] = size(ca{j});
  for i=1:n
     if (j < 7) col = colors(j);
     else col='k'; end;
     
     plot(ca{j}(i,:),col);
     hold on;
  end
end



end
