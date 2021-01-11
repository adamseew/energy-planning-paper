
pow = cswread('');

delta = pow(2,1)-pow(1,1);

for sens=transpose(pow)
   
end

%% model building function

function [A, C] = build_model(omega,r)

    m = 2*r+1;

    Aj = @(omega, j) [0 omega*j ; -omega*j 0];
    A  = zeros(m);
    A(1,1) = 0;
    C = [1];

    for i = 1:r
        A(2*i : 2*i + 1, 2*i : 2*i + 1) = Aj(omega, i); 
        C = [C 1 0];
    end
    
    C = 1/((2*pi)/omega) * C;

end
