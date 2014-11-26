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

function result = noize_dispersion(noized, clear)
  noize = noized - clear;
  result = max(noize) - min(noize);
endfunction

function result = deconvolution(signal, ir)
  figure(); plot2d(ir); xtitle("Impulse Response");
  iri = ir_inversion(ir, 5 * size(ir, 2));
  figure(); plot2d(iri); xtitle("Deconvolution Filter");
  ir_convolution = convol(ir, iri) // should be [1, 0, ..., 0]
  figure(); plot2d(ir_convolution); xtitle("Convolution between filter and reverse filter");
  // Deconvolution is just convolution with reverse-filter, so:
  result = convol(signal, iri);
endfunction

function result = initial_signal(t)
  result = exp(-t);
endfunction

function check_roots(impulse_responce)
  ir_polynom = poly(impulse_responce, 'x', 'coeff');
  ir_roots = roots(ir_polynom);
  ir_roots_modules = (real(ir_roots) .^ 2 + imag(ir_roots) .^ 2) .^ 0.5;
  if find(ir_root_modules < 1) then
    disp("You have roots of impulse response that smaller than 1. Bad news, filter will be unstable");
  end
endfunction

impulse_responce = [0.131, 0.229, 0.268, 0.211, 0.111, 0.039, 0.009, 0.001];
sigma = 0.005; // dispersion of normal noise
delta_t = 0.001; // width between two points, 1 / discretization frequence
t = 0:delta_t:5; // signal points

// Create discrete signal with noize.
//
signal = initial_signal(t)
noize = grand(1, size(t, 2), 'nor', 0, sigma);
noized_signal = signal + noize;
convolved_signal = convol(noized_signal, impulse_responce);

// Deconvolution is just convolution with reverse-filter, so:
deconvolved_signal = deconvolution(convolved_signal, impulse_responce);
deconvolved_signal = deconvolved_signal(1:size(t, 2));
figure(); plot2d(t, deconvolved_signal); xtitle("Deconvolved signal");
printf("Noize dispersion: %f", noize_dispersion(deconvolved_signal, signal));
if (noize_dispersion(deconvolved_signal, signal) < (sigma * 10)) then
  disp("Gain of noize dispersion less then 10: OK");
else
  disp("Gain of noize dispersion less then 10: FAILED");
end
