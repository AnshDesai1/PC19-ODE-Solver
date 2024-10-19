function h=getSteps(ord,t,j)
    for m=1:(ord-1)
        h(m)=t(j-(m-1))-t(j-m); % Take step size between each t(j) and t(j-1)
    end
end