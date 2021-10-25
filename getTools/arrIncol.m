function arrOut = arrIncol(arrIn)
% Return a column array all the time 
% becasue array are so much nicer to print in console
sa = size(arrIn);

if sa(1) == 1
   arrOut = arrIn';
else
    arrOut = arrIn;
end
end