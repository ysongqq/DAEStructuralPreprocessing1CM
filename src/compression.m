function[U,S]=compression(A,p)
sz=size(A);
n=sz(1);
m=sz(2);
U = eye(n);
S=zeros(1,n);

[~,perm] = sort(-p);
inv=perm;
for i=1:n
    inv(perm(i))=i;
end
B=A;


for i=1:n
    B(i,:)=A(perm(i),:);
end

for i=1:n
    piv_j = -1;
        for j=1:m
            if B(i, j) ~= 0
                piv_j =j;
                break
            end
        end
    if piv_j == -1
        S(perm(i))=1;
        continue
    end
    for h = i+1:n
        w = -B(h, piv_j) / B(i, piv_j);
        B(h,:) = B(h,:)+w.*B(i,:);
        U(perm(h),:) = U(perm(h),:)+w.*U(perm(i),:);
    end
end