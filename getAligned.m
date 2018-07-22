function [aligned] = getAligned(bpmH,bpmS,npoints)

[n,m] = size(bpmS);
alignedH2 = zeros(n,npoints);
nal = 0;
for i=1:n
 ix = find(bpmS(i,:) < -1, 1, 'first');
 if (length(ix) > 0)
   nal = nal + 1;
   alignedH2(nal,:) =  bpmH(i,ix:ix+npoints-1);
 end  
end

aligned = alignedH2(1:nal,1:npoints);

end
