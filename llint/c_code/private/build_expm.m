lapacklib = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft','libmwlapack.lib');
blaslib = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft','libmwblas.lib');

mex('-v', '-largeArrayDims', 'expm64v41.c', blaslib, lapacklib)
file = sprintf('expm64v41.%s',mexext);
[status,message,messageId] = movefile(['.',filesep,file],['.',filesep,'..',filesep,file ],'f');
