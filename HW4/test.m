% test to understand fft

t = linspace(0,1000,1001);
x = 2.*cos(t./(2*pi)) + cos(5.*t./(2*pi));

figure
plot(t,x);

f = fft(x);

figure
%semilogy(t,f,'-b');
semilogy(t,abs(real(f)),'-b')

% semilogy(frequency,abs(real(f)),'-b')
% hold on
% semilogy(frequency,abs(imag(f)),'-m');
% xlabel('cycles per day (cpd)');
% xlim([0 125]);

