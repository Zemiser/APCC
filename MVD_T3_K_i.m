function [VectorKStar_i, Flag, zStar] = MVD_T3_K_i(VectorK_i,r,K,N,d,L,alpha,mu,KprimeMax)

Flag = 1;
zMVD = 0;
zStar = (alpha + 1/mu)*K;
VectorK_i_Origin = zeros(1,r);

while sum(abs(VectorK_i - VectorK_i_Origin)) ~= 0
    VectorK_i_Origin = VectorK_i;
    VectorEta_i = d*(VectorK_i + L -1) +1;
    Vector_t_iStar = zeros(1,r);
    for index3 = 1:r
        Vector_t_iStar(index3) = alpha*index3 + index3/mu/(N - VectorEta_i(index3) + 1);
    end
    [zStar,index_j] = max(Vector_t_iStar);
    if VectorK_i(index_j) > 2
        VectorK_iMVD = VectorK_i;
        zMVD = zStar;
        for index_l = 1:r
            if index_l == index_j
                continue;
            end
            Vector_t_iStar_temp = Vector_t_iStar;
            VectorK_i_temp = VectorK_i;
            VectorK_i_temp(index_j) = VectorK_i(index_j) - 1;
            Vector_t_iStar_temp(index_j) = alpha*index_j + index_j/mu/(N - (d*(VectorK_i_temp(index_j) + L -1) +1) + 1);
            if VectorK_i_temp(index_l) <= KprimeMax
                VectorK_i_temp(index_l) = VectorK_i(index_l) + 1;
                Vector_t_iStar_temp(index_l) = alpha*index_l + index_l/mu/(N - (d*(VectorK_i_temp(index_l) + L -1) +1) + 1);
                [zStar_temp,index_j_temp] = max(Vector_t_iStar_temp);
                if zStar_temp < zMVD 
                    VectorK_iMVD = VectorK_i_temp;
                    zMVD = zStar_temp;
                end
                if zStar_temp == zMVD && index_l ~= index_j_temp
                    VectorK_iMVD = VectorK_i_temp;
                    zMVD = zStar_temp;
                end
            end
        end
        VectorK_i = VectorK_iMVD;
    else
        VectorKStar_i = VectorK_i;
        break;
    end
end
VectorKStar_i = VectorK_i;

end