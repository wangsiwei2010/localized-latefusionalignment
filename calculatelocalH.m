function [localHmatrix] = calculatelocalH(H,A)
[num,classnumber]= size(H);
localHmatrix = zeros(num,classnumber);
for i=1:num
    localHmatrix(i,:) = H(i,:)*A(i,2);
end

localHmatrix = (1/num) * localHmatrix;

end