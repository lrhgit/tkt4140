%program gtest
nrpm = 6000
for v = [13.7 21.6 29.9 38.4 46.9 55.2 63.1 71.9 80.2 88.1]
    [cd , cl] = cdcldata(v,nrpm);
    fprintf( '%7.2f %10.4f %10.4f \n',v,cd,cl);
end
