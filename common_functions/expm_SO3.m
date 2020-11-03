function R = expm_SO3(r)
theta=norm(r);

y=sinx_over_x(theta);
y2=sinx_over_x(theta/2);

R=eye(3)+y*hat(r)+1/2*y2^2*hat(r)^2;
end