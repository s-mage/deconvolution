// TODO:
// h^-1 (n) = (1 / h(0)) sum_{k = 0}^{n - 1} h^-1 (k) h(n - k)

function initial_signal(t)
  exp(-t);
endfunction

impulse_responce = [0.131, 0.229, 0.268, 0.211, 0.111, 0.039, 0.009, 0.001];
sigma = 0.005; // dispersion of normal noise
delta_t = 0.1; // width between two points, 1 / discretization frequence
t = 0:delta_t:5; // signal points

ir_polynom = poly(impulse_responce, 'x', 'coeff');
ir_roots = roots(ir_polynom);
ir_roots_modules = (real(ir_roots) .^ 2 + imag(ir_roots) .^ 2) .^ 0.5;

function ir_inversion(ir, n)
  result = zeros(1, n);
  result(1) = 1 / ir(1);
  for k = 1:(n - 1)
    ir_limited = ir(1:k);
    ir_limited_reversed = ir_limited(k:-1:1)'; // this is how to reverse vector in scilab, no matter.
    result(k + 1) = 1 / ir(1) * (result(1:k) * ir_limited_reversed);
  end
endfunction
