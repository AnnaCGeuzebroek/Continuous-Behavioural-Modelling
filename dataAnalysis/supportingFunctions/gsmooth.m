function y = gsmooth(x,fs,sd,sdextra)
% Gaussian smoothing
% x  = data (samples)
% y  = smoothed data
% fs = sample frequency data (samples/s = Hz)
% sd = standard deviation gauss (in seconds)
% sdextra = number of sd's added to data edges (5 = default)
%
% BeerendW, 

if nargin<4
    sdextra = 5;
end

x           = x(:)';                % Make sure data are in a row vector
nextra      = round(sdextra*sd*fs); % Extra sample points?
nx          = length(x);            % Length(x)

% Number of fft-points + 2 times extra sd's
nfft        = nx+2*nextra;
% determine frequency for each point
if rem(nfft,2)                                  % nfft odd
    frq     = [0:(nfft-1)/2 (nfft-1)/2:-1:1];
else
    frq     = [0:nfft/2 nfft/2-1:-1:1];
end
% Add extra sd
x = [x(nextra+1:-1:2) x x(end:-1:end-nextra+1)]; % mirror edges to go around boundary effects
% Determine gaussian
g = exp(-0.5 * ((frq - 0)./sd*fs).^2) ./ (sqrt(2*pi) .* sd*fs);
g = g./sum(g);
% filter the signal with fast fourier 
y = real(ifft(fft(x).*fft(g)));
y = y(nextra+1:end-nextra);
