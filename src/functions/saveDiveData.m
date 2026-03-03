function saveDiveData(filePath, Dives)
%SAVEDIVEDATA Summary of this function goes here

tic
fprintf('Caching processed dives to %s...\n', filePath);
save(filePath, 'Dives','-v7.3');
toc
fprintf('Output saved.\n');

end