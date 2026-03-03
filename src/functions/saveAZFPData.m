function saveAZFPData(savePath, Output)
%SAVEAZFPDATA Caches output data in file

tic
fprintf('Caching processed output to %s...\n', savePath);
save(savePath, 'Output', '-v7.3');
fprintf('Output saved and will not be re-generated if not deleted.\n');
toc
end
