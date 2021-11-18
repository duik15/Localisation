function mkdir2(folder)

if ~isfolder(folder) 
    disp(['Creating output folder: ' folder]);
    mkdir(folder); 
end
end