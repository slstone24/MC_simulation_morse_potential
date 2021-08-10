function E = Energymorse(r)

D = 1; % well depth
a = 1; % well width
re = 1; % eq bond distance

% compute energy over grid: sum of Natom pairwise interactions
Natom = size(r,1);
E = 0;
    for j=1:Natom
        for k=j+1:Natom
            rjk = norm(r(j,:)-r(k,:));
            E = E + D*(1-exp(-a*(rjk-re)))^2;
        end
    end
    
end