function [u] = riemann(ul, ur, s)

if(ul > ur) %Shockwave
    if(0 < s(ul, ur)) %Right shockwave
        u = ul;
    else %Left shockwave
        u = ur;
    end
else %Rarefaction wave
    if(ul > 0) %Right expansion
        u = ul;
    elseif(ur < 0) %Left expansion
        u = ur;
    else % Centered expansion
        u = 0;
    end
end

end