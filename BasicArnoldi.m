function [v,H] = BasicArnoldi(A,v1,m)
% Computes the basic Arnoldi algorithm for a matrix A and a vector v1
% The parameter m specifies the maximum order of the computed Krylov subspace
% Stores the Krylov subspace basis vectors in v and the Hessenberg matrix in H

    v(:,1) = v1/norm(v1);
    for j=1:m
        w=A*v(:,j);
        for i=1:j
            H(i,j)=w'*v(:,i);
            w = w - H(i,j)' * v(:,i);
        end
        H(j+1,j)=norm(w);
        if H(j+1,j)<1e-5 %Stop at 0
            H = H(1:j,1:j);
            break
        end
        v(:,j+1) = w/H(j+1,j);
    end
end
