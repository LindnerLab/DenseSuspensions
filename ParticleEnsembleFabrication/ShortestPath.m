function nodesSorted = ShortestPath(nodes)
nodesSorted = [];
Nnodes = length(nodes);
p = [0 ,0];
for i = 1:Nnodes
    d = sqrt(sum((nodes-p).^2,2));
    [~, idx] = min(d);
    p = nodes(idx,1:2);
    nodesSorted = [nodesSorted; p];
    nodes(idx,:) = [];
end
end