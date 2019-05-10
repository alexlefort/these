function [out] = interpolation(p_n,p_xt,p_yt,p_x)

y = 0;
if( p_n == 0 )
    y=0.0;
elseif (p_n == 1)
    y=p_yt(1);
elseif (isnan(p_x))
    y=p_x;
elseif (p_x <= p_xt(1))
    y=p_yt(1);
elseif (p_x >= p_xt(p_n))
    y=p_yt(p_n);
else
    ii = 1;
    il = 1;
    ih = p_n;
    a = 0;
    while ((ii<=p_n) && (ih - il > 1))
        ii = ii + 1;
        im=floor((il+ih)/2);
        if(p_xt(im) >= p_x) ih=im;
        end
        if(p_xt(im) <= p_x) il=im;
        end
    end
    if (il == ih)
        a =0.0;
    else
        a = (p_xt(ih)-p_x)/(p_xt(ih)-p_xt(il));
    end
    y = a*p_yt(il) + (1.0-a)*p_yt(ih);
end
out = y;
