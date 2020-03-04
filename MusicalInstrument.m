% Musical instrument



for clock = 1:clockmax

    t = clock*dt;
    j = 2:J-1;
    V(j) = V(j) + (dt/dx^2) * (T/M) * (H(j+1)-2*H(j)+H(j-1)); % vector operation
    H(j) = H(j) + dt * V(j);
end
