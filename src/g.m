% error function sigmoid    

function out = g(v,v0,varsigma)

    out = 0.5*erf((v - v0) / (sqrt(2)*varsigma)) + 0.5;
%       out = 1 ./ (1 + exp(varsigma*(-v+v0)));
    
end