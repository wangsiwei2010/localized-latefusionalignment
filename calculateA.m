function [A] = calculateA(NS,neibournum)
[~,samplenumber] = size(NS);
A = zeros(samplenumber,2);
A(:,1)=[1:1:samplenumber]';
for i=1:neibournum
    for j=1:samplenumber
        a = NS(i,j);
        A(a,2) = A(a,2)+1;
    end
end

end

