function aOut = shiftAngle90(aIn)
aOut = aIn -90;
aOut(aOut<0) = 360  + aOut(aOut<0)
end