% Draw 3-D Antenna
close all
N = 28;
Amp = dx(1:N);
Phi = dx(N+1:2*N);
Pos = dx(2*N+1:end);
subplot(1,2,1);
p09_antenna(x,'draw');
subplot(1,2,2);
hold on
for i = 1:N
    plot3([Pos(i),Pos(i)],[0,Amp(i)*sin(Phi(i))],[0,Amp(i)*cos(Phi(i))]);
end