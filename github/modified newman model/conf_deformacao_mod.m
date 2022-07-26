function Li_max = conf_deformacao_mod(d_max,xi_max,wi_max,Li_max_ret,n_zone_max,xi,Li_max,n_zone)

for i=1:n_zone

        if xi(i)<(d_max)
            for j=n_zone_max:-1:1
                if (xi_max(j)-wi_max(j))<xi(i) && xi(i)<(xi_max(j)+wi_max(j))
                    if Li_max(i) < Li_max_ret(j)
                    Li_max(i) = Li_max_ret(j);
                    end
                    break;
                elseif (xi_max(j)-wi_max(j))==xi(i)
                    if j~=n_zone_max %esse caso de j=10 e redundante pq e impossivel de acontecer, mas mantive para manter o raciocinio 
                        if Li_max(i)<(Li_max_ret(j)+Li_max_ret(j+1))/2
                        Li_max(i) = (Li_max_ret(j)+Li_max_ret(j+1))/2; %aq e j+1 pq e o elemento de tras q compartilha a borda de elemento
                        end
                    else
                        if Li_max(i) < Li_max_ret(j)
                        Li_max(i) = Li_max_ret(j);
                        end
                    end
                    break;
                elseif xi(i)==(xi_max(j)+wi_max(j))
                    if j~=1
                        if Li_max(i)<(Li_max_ret(j)+Li_max_ret(j-1))/2
                        Li_max(i) = (Li_max_ret(j)+Li_max_ret(j-1))/2; %aq e j-1 pq e o elemento a frente q compartilha a borda de elemento
                        end
                    else
                        if Li_max(i) < Li_max_ret(j)
                        Li_max(i) = Li_max_ret(j);
                        end
                    end
                    break;
                end
            end
        elseif xi(i)==(d_max)
            Li_max(i) = (Li_max(i)+Li_max_ret(i))/2;
        end
        
        %% mesmo codigo mas sem a condicao de manter a maxima deformacao atual, caso ela seja maior
        
%          if xi(i)<(d_max)
%             for j=n_zone_max:-1:1
%                 if (xi_max(j)-wi_max(j))<xi(i) && xi(i)<(xi_max(j)+wi_max(j))
%                     %if Li_max(i) < Li_max_ret(j)
%                     Li_max(i) = Li_max_ret(j);
%                     %end
%                     break;
%                 elseif (xi_max(j)-wi_max(j))==xi(i)
%                     if j~=n_zone_max %esse caso de j=10 e redundando pq e impossivel de acontecer, mas mantive para manter o raciocinio 
%                         %if Li_max(i)<(Li_max_ret(j)+Li_max_ret(j+1))/2
%                         Li_max(i) = (Li_max_ret(j)+Li_max_ret(j+1))/2; %aq e j+1 pq e o elemento de tras q compartilha a borda de elemento
%                         %end
%                     else
%                         %if Li_max(i) < Li_max_ret(j)
%                         Li_max(i) = Li_max_ret(j);
%                         %end
%                     end
%                     break;
%                 elseif xi(i)==(xi_max(j)+wi_max(j))
%                     if j~=1
%                         %if Li_max(i)<(Li_max_ret(j)+Li_max_ret(j-1))/2
%                         Li_max(i) = (Li_max_ret(j)+Li_max_ret(j-1))/2; %aq e j-1 pq e o elemento a frente q compartilha a borda de elemento
%                         %end
%                     else
%                         %if Li_max(i) < Li_max_ret(j)
%                         Li_max(i) = Li_max_ret(j);
%                         %end
%                     end
%                     break;
%                 end
%             end
%         elseif xi(i)==(d_max)
%             Li_max(i) = (Li_max(i)+Li_max_ret(i))/2;
%         end
        
% se a coordenada do elemento nao estiver nem dentro da zona max e se sua
% coordenada nao coincidir com a borda do elemento mais extremo da zona
% max, entao o elemento esta considerado fora da zona plastica maior e
% mantem seu valor ja calculado de Li

    end
end