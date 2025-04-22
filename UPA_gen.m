function a=UPA_gen(N)
theta=2*pi*rand(1);
phi=pi*rand(1)-pi/2;
N_x=floor(sqrt(N));
a=[];
for n=1:N
    phas=pi*(floor(n/N_x)*sin(theta)*sin(phi)+(n-floor(n/N_x)*N_x)*sin(phi)*cos(theta));
    temp=exp(j*phas);
    a=[temp;a];
end