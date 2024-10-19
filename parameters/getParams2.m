function [pbeta,cbeta,err]=getParams2(h,H)
    pbeta=[
    -H*H/(2*h(1))
    H+H^2/(2*h(1))
    ];
    
    cbeta=[
    0.5*H
    0.5*H
    ];
    
    err=[
    1/12*H^2*(2*H+3*h(1))
    -1/12*H^3
    ];
end