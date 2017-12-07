a=randn(100,1000)+ cos(2*pi/10*(1:100)')*ones(1,1000);
b=randn(100,1000) + sin(2*pi/10*(1:100)')*ones(1,1000);

fa=fft(a);
fb=fft(b);
fab=conj(fa).*fb;
faa=conj(fa).*fa;
fbb=conj(fb).*fb;

cab=abs(mean(fab,2)) ./sqrt(abs(mean(faa,2)) .* abs(mean(fbb,2)));

m=10;
clear phase_c
for i=1:1000/m
phase_c(:,i)=atan2(-imag(mean(fab(:,(i-1)*m+1:i*m),2)),...
real(mean(fab(:,(i-1)*m+1:i*m),2)));
end

nd=m;
delta_phase = asin(tinv(.95,2*nd)*...
sqrt((1-abs(cab).^2)./(abs(cab).^2*sqrt(2*nd))));
delta_phase2 = sqrt((1-cab.^2)./(abs(cab).^2*2*nd));
delta_phase3 = asin(tinv(.975,2*nd-2)*(1 ./cab.^2-1)/(2*nd-2));
% compare results
[delta_phase(11) delta_phase2(11) delta_phase3(11)] % we like delta_phase2 best
std(phase_c(11,:))