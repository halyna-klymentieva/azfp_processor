function saveDiveData(filePath, Dives, bottomDepth)
%SAVEDIVEDATA Summary of this function goes here

tic
fprintf('Caching processed dives to %s...\n', filePath);
save(filePath, 'Dives', 'bottomDepth', '-v7.3');
toc
fprintf('Output saved.\n');

end