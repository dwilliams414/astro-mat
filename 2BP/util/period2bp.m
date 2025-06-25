function period = period2bp(a, GM)
%PERIOD2BP: Compute orbit period from semimajor axis (a) and mass parameter
arguments
    a (1, :) double
    GM (1, :) double
end
    if length(GM) == 1
        GM = GM.*ones(size(a));
    end
    
    period = 2*pi*sqrt((a.^3)./GM);
end

