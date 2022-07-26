function Leq = defor_eq(crack,Li_max,wi)

n = length(Li_max);
bar_length = 2*wi;

if crack <= bar_length(n)
    Leq = Li_max(n);
else
    for i=n-1:-1:1
        if crack <= sum(bar_length(i:n))
            somatorio1 = 0;
            somatorio2 = 0;
                for j=n:-1:i
                    somatorio1 = somatorio1 + Li_max(j)*wi(j);
                    somatorio2 = somatorio2 + wi(j);
                end
            Leq = somatorio1/somatorio2;
            break;
        end
    end
end

end