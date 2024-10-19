function nextH=estimateNextH(err,ord,u1,u2,etarg,h)
    d=norm((u2-u1)/(err(1)-err(2))); % Compute estimate for derivative
    % Acquire variables for Newton's Method
    [E,dE]=choosefunc(ord,d,h,etarg); % Get function for error and derivative
    nHc=[12 24 720/19 160/3 302400/4315 4233600/48125 101606400/950684 1036800/8138 914457600/7217406]; % nextH coefficients
    nextH=(nHc(ord-1)*etarg/d)^(1/ord); % initial nextH estimate
    lambda=1e-4; % Error Ratio
    epsilon=1e-14; % Min error
    E0=min(abs(E(nextH)),1e-4); % Initial error
    E1=1; % Max error
    its=0; % Starting iteration
    maxits=20; % Maximum iterations
    while its<maxits && E1>lambda*E0 && E1>epsilon % Newton's Method to estimate nexth 
        nextH=nextH-E(nextH)/dE(nextH);
        E1=abs(E(nextH));
        its=its+1;
    end
end