// TODO:
// h^-1 (n) = (1 / h(0)) sum_{k = 0}^{n - 1} h^-1 (k) h(n - k)

function result = initial_signal(t)
  result = exp(-t);
endfunction

impulse_responce = [0.131, 0.229, 0.268, 0.211, 0.111, 0.039, 0.009, 0.001];
sigma = 0.005; // dispersion of normal noise
delta_t = 0.1; // width between two points, 1 / discretization frequence
t = 0:delta_t:5; // signal points

ir_polynom = poly(impulse_responce, 'x', 'coeff');
ir_roots = roots(ir_polynom);
ir_roots_modules = (real(ir_roots) .^ 2 + imag(ir_roots) .^ 2) .^ 0.5;

function result = ir_inversion(ir, n)
  result = zeros(1, n);
  result(1) = 1 / ir(1);
  ir = ir(2:$);
  N = size(ir, 2);
  for k = 1:(n - 1)
    if(k < N) then h = ir(k:-1:1); else h = ir(:,$:-1:1); end
    if(k < N) then hi = result(1:k); else hi = result(k - N + 1:k); end
    result(k + 1) = - result(1) * (hi * h');
  end
endfunction

ir = impulse_responce;
figure();
plot2d(ir);
iri = ir_inversion(ir, 3 * size(ir, 2))
figure();
plot2d(iri);
ir_convolution = convol(ir, iri) // should be [1, 0, ..., 0]
figure();
plot2d(ir_convolution);

// Create discrete signal with noize.
//
signal = initial_signal(t)
noize = grand(1, size(t, 2), 'nor', 0, sigma);
noized_signal = signal + noize;

convolved_signal = convol(noized_signal, impulse_responce);

// Deconvolution is just convolution with reverse-filter, so:
deconvolved_signal = convol(convolved_signal, iri);

figure();
plot2d(t, deconvolved_signal(1:51));
