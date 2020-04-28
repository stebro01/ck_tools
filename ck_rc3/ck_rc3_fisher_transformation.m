function Z = ck_rc3_fisher_transformation(R)
% Fishers Z-Transformation
% für Matrizen beliebiger Dimension
% z = 0.5 * log[(1 + r)/(1-r)]
%

Z = zeros(size(R));
dim = size(size(R),2);

switch dim
    case 1
        for i=1:size(R,1)
            y = R(i);
            Z(i)=0.5 * log((1+y)/(1-y));
        end
    case 2
        for i=1:size(R,1)
            for j=1:size(R,2)
                y = R(i,j);
                Z(i,j)=0.5 * log((1+y)/(1-y));
            end
        end

    case 3
        for i=1:size(R,1)
            for j=1:size(R,2)
                for k=1:size(R,3)
                    y = R(i,j,k);
                    Z(i,j,k)=0.5 * log((1+y)/(1-y));
                end
            end
        end

    case 4
        for i=1:size(R,1)
            for j=1:size(R,2)
                for k=1:size(R,3)
                    for l=1:size(R,4)
                        y = R(i,j,k,l);
                        Z(i,j,k,l)=0.5 * log((1+y)/(1-y));
                    end
                end
            end
        end
end
end