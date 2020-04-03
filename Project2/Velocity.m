% Calculate velocity from distance to front car
function [res] = Velocity(d, dmin, dmax, vmax)

if(d < dmin)
    res=0;
elseif(d < dmax)
    res = vmax * log(d/dmin) / log(dmax/dmin);
else
    res = vmax;
end

end
