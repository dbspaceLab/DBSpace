function tri = triangle(t)
%TRIANGLE Gives full triangle wave (upslope and downslope)
%   phase works just like a sine wave
    tri = 2*(abs(sawtooth(t))-0.5);
end

