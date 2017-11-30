function [ z ] = Broaden( x,y,sigma )
%Broaden Broadens a spectral line with a gaussian
%    z  = Broaden( x,y,sigma )
%       z = output broadened line
%       x = input energy grid
%       y = input spectral line
%       sigma = Gaussian sigma

L = max(x)-min(x);
x2 = linspace(-L/2,L/2,numel(x));
cony = exp(-x2.^2/(2*sigma^2));
cony = cony/trapz(x2,cony);

Norm = trapz(x,y);
z = conv(y,cony,'same');
z = z/trapz(x,z)*Norm;

end

