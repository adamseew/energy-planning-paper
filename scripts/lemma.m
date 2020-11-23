

% using this script to temporarily proof the lemma

a0 = 2000;
%a0 = qt(1);
a1 = 200;
%a1 = qt(2);

b1 = 200;
%b1 = qt(3);

a2 = 200;
b2 = 200;

a3 = 200;
b3 = 200;

T = 100;

o = 2*pi/T;

h = @(k) a0/T + 2*(a1*cos(o*k) + b1*sin(o*k))/T + ...
                2*(a2*cos(o*2*k) + b2*sin(o*2*k))/T + ...
                2*(a3*cos(o*3*k) + b3*sin(o*3*k))/T;
              
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

Aj = @(j) [ 0 o*j ; -o*j 0 ];
A  = zeros(7);
A(1,1) = 0;
C = 1;

for i = 1:3
    A(2 * i : 2 * i + 1, 2 * i : 2 * i + 1) = Aj(i); 
    C = [C 1 0];
end

q0 = [a0; 2*a1; 2*b1; 2*a2; 2*b2; 2*a3; 2*b3]; 
C = 1/T * C;

yy = [];
    
q = [];

%A(2,3) = 0;
%A(3,2) = 0;
%A(3,3) = A(3,2);

for k = 1:.01:x(end)
    q1 = (eye(size(q0,1)) + A*.01)*q0;
    
    y1 = C * q1;
    
    yy = [yy y1];
    q = [q q1];

    q0 = q1;
end

figure;
plot(1:.01:x(end), yy);

yyy = c(x);

figure;
plot(x, yyy);

