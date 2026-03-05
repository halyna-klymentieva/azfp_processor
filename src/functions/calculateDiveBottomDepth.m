function bottomDepth = calculateDiveBottomDepth(azfpData, Dives)
%CALCULATEDIVEBOTTOMDEPTH Calculate bottom depth and end dives
nDives = length(Dives);
endDive = zeros(nDives, 1);
bottomDepth = zeros(nDives, 1);

for i = 1:nDives
    endDive(i) = Dives(i).Index(2) - 1;
    bottomDepth(i) = azfpData(1).Depth(endDive(i));
end

end