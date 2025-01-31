%% Euler's formula
% formula: M * e^(i*k) = M (cos(k) + i*sin(k))
% i : the imagery operator. square root of i is -1
% k : angle with respect to the positive real axis (phase)
% M : length or distance from the origin (magnitude)
% this formula describe a vector in a complex (polar) plane.
% will be very useful later on interpreting fourier coefficients.
M = 10;
k = pi/6;

meik = M*exp(1i*k);

%plot on polar plane
figure, clf
subplot(121)
polar([0 k],[0 M],'r'), hold on
polar(k, M,'ro')
title('Polar plane')

%plot on complex plane
subplot(122), hold on
plot(meik,'ro')
plot(real(meik),imag(meik),'gs')
axis([-1 1 -1 1]*abs(meik))
axis square
xlabel('Real'), ylabel('Imag')
grid on
title('Cartesian (rectangular) plane')