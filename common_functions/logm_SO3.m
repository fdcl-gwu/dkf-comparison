function r = logm_SO3(R)

[V lam]=eig(R);
eps=1e-6;
min_del_lam_1=1;
for i=1:3
    if norm(imag(V(:,i))) < eps
        if (lam(i,i)^2-1) < min_del_lam_1
            min_del_lam_1=lam(i,i)^2-1;
            i_min=i;
        end
    end
end
v=real(V(:,i_min));

cos_theta=(trace(R)-1)/2;
if cos_theta > 1.0
    cos_theta=1;
elseif cos_theta < -1
    cos_theta=-1;
end
theta=real(acos(cos_theta));
R_new=expm_SO3(theta*v);

if norm(R-R_new) > norm(R-R_new')
    v=-v;
end

r=v*theta;
end