

% using this script to temporarily proof the lemma

a0 = 10;
a0 = qt(1);
a1 = 1.43;
a1 = qt(2);
a2 = 1.43;
a2 = qt(3);
a3 = 1.43;
a3 = qt(4);
b1 = 1.43;
b1 = qt(5);
b2 = 1.43;
b2 = qt(6);
b3 = 1.43;
b3 = qt(7);
T = 100;

o = 2*pi/T;

h = @(k) a0/T + 2*(a1*cosd(o*k) + b1*sind(o*k))/T + ...
                2*(a2*cosd(o*2*k) + b2*sind(o*2*k))/T + ...
                2*(a3*cosd(o*3*k) + b3*sind(o*3*k))/T;
              
c = @(k) a0/T + (exp(1i*o*k)*(a1 - 1i*b1))/T + ...
                (exp(-1i*o*k)*(a1 + 1i*b1))/T + ...
                (exp(1i*o*2*k)*(a2 - 1i*b2))/T + ...
                (exp(-1i*o*2*k)*(a2 + 1i*b2))/T + ...
                (exp(1i*o*3*k)*(a3 - 1i*b3))/T + ...
                (exp(-1i*o*3*k)*(a3 + 1i*b3))/T;
              
x = 1:1000;
y = h(x);
              
figure;
plot(x, y);

Aj = @(j) [ 1 o*j ; -o*j 1 ];
A  = zeros(7);
A(1,1) = 1;
C = 1;

for i = 1:3
    A(2 * i : 2 * i + 1, 2 * i : 2 * i + 1) = Aj(i); 
    C = [C 1 1];
end

q0 = [a0; a1; b1; a2; b2; a3; b3];    
C = 1/T * C;

yy = [];
    
q = [];
    
for k = x
    q1 = A * q0;
    y1 = C * q1;
        
    yy = [yy y1];
    q = [q q1];
        
    q0 = q1;    
end

figure;
plot(x, yy);

yyy = c(x);


figure;
plot(x, yyy);

