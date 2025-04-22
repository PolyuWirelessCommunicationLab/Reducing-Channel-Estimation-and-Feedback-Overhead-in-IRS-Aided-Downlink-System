function a=ULA_gen(M)
theta=2*pi*rand(1);
a=[];
for m=1:M
    phas=pi*(m-1)*sin(theta);
    temp=exp(j*phas);
    a=[temp;a];
end