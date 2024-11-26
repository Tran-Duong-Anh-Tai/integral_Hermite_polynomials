% Creating dictionary Wijkl, Yijkl, Uijklmn
Wijkl   = dictionary();
Yijkl   = dictionary();
Uijklmn = dictionary();

% Input value Morb
Morb = 6;

% Printing W_{0,0,0,0}
Wijkl("0000") = 1 / sqrt(2 * pi);

% Calculating Wijkl
for i = 1:Morb
    for j = 0:i
        for k = 0:j
            for l = 0:k
                if mod(i + j + k + l, 2) == 0
                   % Creating index
                   index = strcat(string(i), string(j), string(k), string(l));

                   % Calculating the value of the integral using a recursive formula
                   Wijkl(index) = 0.5*(-sqrt((i - 1) / i) * calc_W(Wijkl, [i - 2, j, k, l]) + ...   
                                        sqrt(j / i) * calc_W(Wijkl, [i - 1, j - 1, k, l])   + ...
                                        sqrt(k / i) * calc_W(Wijkl, [i - 1, j, k - 1, l])   + ...
                                        sqrt(l / i) * calc_W(Wijkl, [i - 1, j, k, l - 1]));
                end
            end
        end
    end
end
disp(Wijkl);

% Calculating Yijkl
for i = 0:Morb
    for j = 0:i
        for k = 0:j
            for l = 0:k
                if mod(i + j + k + l, 2) == 0 && ( (i ~= j) && (k ~= l) )
                   % Creating index
                   index = strcat(string(i),string(j),string(k),string(l));

                   % Calculating the value of the integral using a recursive formula
                   Yijkl(index) = 2*(sqrt(i * k) * calc_W(Wijkl, [i - 1, j, k - 1, l]) - ...   
                                     sqrt(i * l) * calc_W(Wijkl, [i - 1, j, k, l - 1]) - ...
                                     sqrt(j * k) * calc_W(Wijkl, [i, j - 1, k - 1, l]) + ...
                                     sqrt(j * l) * calc_W(Wijkl, [i, j - 1, k, l - 1]));
                  
                end
            end
        end
    end
end
disp(Yijkl);  

% Printing U_{0,0,0,0,0,0}
Uijklmn("000000") = 1 / (pi*sqrt(3));

% Calculating Uijklmn
for i = 1:Morb
    for j = 0:i
        for k = 0:j
            for l = 0:k
                for m = 0:l
                    for n = 0:m
                        if mod(i + j + k + l + m + n, 2) == 0                     
                            % Creating index
                            index = strcat(string(i),string(j),string(k),string(l),string(m),string(n));

                            %  Calculating the value of the integral using a recursive formula
                            Uijklmn(index) = (1/3)*(-2*sqrt((i - 1) / i) * calc_U(Uijklmn, [i - 2, j, k, l, m, n]) + ...   
                                                       sqrt(j / i) * calc_U(Uijklmn, [i - 1, j - 1, k, l, m, n])   + ...
                                                       sqrt(k / i) * calc_U(Uijklmn, [i - 1, j, k - 1, l, m, n])   + ...
                                                       sqrt(l / i) * calc_U(Uijklmn, [i - 1, j, k, l - 1, m ,n])   + ...
                                                       sqrt(m / i) * calc_U(Uijklmn, [i - 1, j, k, l, m - 1 ,n])   + ...
                                                       sqrt(n / i) * calc_U(Uijklmn, [i - 1, j, k, l, m, n - 1]) );
                        end
                    end
                end
            end
        end
    end
end
disp(Uijklmn)





% Function to lookup Wijkl
function result = calc_W(Wijkl,index)
    if max(index) > 0 && min(index) >= 0
        index  = sort(index,'descend');
        result = Wijkl(strcat(string(index(1)),string(index(2)),string(index(3)),string(index(4)))); 
    elseif all(index == 0)
         result = 1 / sqrt(2 * pi);
    else 
         result = 0; 
    end
end

% Function to lookup Uijlkmn
function result = calc_U(Uijklmn,index)
    if max(index) > 0 && min(index) >= 0
        index  = sort(index,'descend');
        result = Uijklmn(strcat(string(index(1)),string(index(2)),string(index(3)),string(index(4)),string(index(5)),string(index(6)))); 
    elseif all(index == 0)
        result = 1 / (pi*sqrt(3));
    else
        result = 0; 
    end
end
