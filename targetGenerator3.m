function v = targetGenerator3(Ps, fdg, t, ri, lambda, tau_p, tau, alpha)
c = 3e8;

l1 = sqrt(Ps) * exp(-1j*2*pi*fdg*t) .* exp(-1j*4*pi*ri/lambda);
l2 = tau_p - abs(tau-2*ri/c);
l3 = sinc(alpha*(tau-2*ri/c) .* (tau_p-abs(tau-2*ri/c)));
l4 = rectpuls( (tau-2*ri/c) / (2*tau_p) );

v = l1.*l2.*l3.*l4;

end