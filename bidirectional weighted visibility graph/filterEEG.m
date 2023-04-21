function y = filterEEG(x,F3dB1,F3dB2)
persistent Hd;

if isempty(Hd)
    
    % The following code was used to design the filter coefficients:
    
    N     = 10;    % Order
    Fs    = 173.61;  % Sampling Frequency
    
    h = fdesign.bandpass('n,f3db1,f3db2', N, F3dB1, F3dB2, Fs);
    
    Hd = design(h, 'butter', ...
        'SystemObject', true);
    

end

s = double(x);
y= step(Hd,s);


