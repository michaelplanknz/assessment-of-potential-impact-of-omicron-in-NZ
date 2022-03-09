function r = calcGRfromReff(Reff, genA, genB)

par.infectKernScale = genA;
par.infectKernShape = genB;

g = wblstat(par.infectKernScale, par.infectKernShape);
a = (1e-6:0.01:50);     % start at slightly >0 to allow for PDFs that are singular at the origin
w = wblpdf(a, par.infectKernScale, par.infectKernShape);
w = w/trapz(a, w);
r = zeros(size(Reff));
myFn = @(x, Re)(Re - 1./trapz(a, w.*exp(-x.*a), 2 ));
xa = (-1:0.01:1)';
for ii = 1:length(Reff)
    r0 = (Reff(ii)-1)/g;    %exponential
    r(ii) = fzero(@(x)myFn(x, Reff(ii)), r0);
end
