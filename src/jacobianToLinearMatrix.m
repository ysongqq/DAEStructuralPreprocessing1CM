function [J0,Jb,Jc,m,n]=jacobianToLinearMatrix(J)
n=size(J,1);
J0=zeros(n,n);
numbers='1234567890.';
syms D
dict=containers.Map;
for i=1:size(J,1)
    %disp(i)
    for j=1:size(J,2)
        f=children(expand(J(i,j)+D,'ArithmeticOnly',true));
        for k=1:numel(f)
            if(char(f(k))=='D')
                continue
            end
            c="";
            isminus=false;
            s=char(f(k));
            if(s=='0')
                break
            end
            for l=1:numel(s)
                if(l==1 && s(1)=='-')
                    isminus=true;
                    continue
                end
                if any(s(l)==numbers)
                    c=append(c,s(l));
                else
                    break
                end
            end
            if(c=="")
                d=1;
            else
                d=double(c);
            end
            if(isminus)
                d=-d;
            end
            g=string(f(k)/d);
            if(g=="1")
                J0(i,j)=J0(i,j)+d;
            elseif(isKey(dict,g))
                mat=dict(g);
                mat(i,j)=mat(i,j)+d;
                dict(g)=mat;
            else
                dict(g)=zeros(n,n);
                mat=dict(g);
                mat(i,j)=mat(i,j)+d;
                dict(g)=mat;
            end 
        end
    end
end
values=dict.values;
keys=dict.keys;
m=0;
syms z1 z2 z3 z4 z5 z6 z7 z8 z9 z10 z11 z12 z13 z14 z15 z16 z17 z18 z19 z20 z21 z22 z23 z24 z25 z26 z27 z28 z29 z30 z31 z32 z33 z34 z35 z36 z37 z38 z39 z40 z41 z42 z43 z44 z45 z46 z47 z48 z49 z50 z51 z52 z53 z54 z55 z56 z57 z58 z59 z60 z61 z62 z63 z64 z65 z66 z67 z68 z69 z70

Jb=[];
Jc=[];
cnt=0;

for i=1:length(values)% rank decomposition
    M=values{i};
    rk=rank(M);
    m=m+rk;
        j=0;
        for k=1:n
            rw=M(:,k);
            if(any(rw~=0))
                j=j+1;
                Jb=[Jb,rw];
            else
                continue
            end
            for l=1:n
                if(rw(l)~=0)
                    cl=M(l,:)/rw(l);
                    Jc=[Jc;cl];
                    M=M-rw*cl;
                    break
                end
            end
        end
end
disp("total rank "+m+", number of coefficient matrices "+length(values))
if(m-length(values)~=0)
    error('The linear symbolic matrix is not a 1CM-matrix.');
end
% preprocessing
mat=zeros(m,n*n);
for i=1:m
    Ji=Jb(:,i)*Jc(i,:);
    mat(i,:)=Ji(:).';
end
[X,Y]=compression(mat,1:m);
Jb=Jb(:,find(Y==0));
Jc=Jc(find(Y==0),:);
m=size(Jc,1);
