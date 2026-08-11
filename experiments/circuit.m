function [elapsed]=circuit()

% CIRCUIT    Apply the proposed method to the transistor amplifier DAE
%            and return the elapsed time in seconds.
%
%   Run "setup" once before calling this function.

[eqs,vars]=problem.transistorAmplifier();

syms s
rowtrans=true;

disp("start")
tic
while(true)
    [J,p,q]=systemJacobian(eqs,vars);
    [J0,Jb,Jc,m,n]=jacobianToLinearMatrix(J);
    % build the sparse representation from the system Jacobian (1CM-matrix)
    Q1=cat(2,cat(2,eye(m),zeros(m,m)),Jc);
    Q2=cat(2,cat(2,zeros(n,m),Jb),J0);
    Q=cat(1,Q1,Q2);
    T=cat(2,cat(2,eye(m),eye(m)),zeros(m,n));
    [rk,dual] = LMMatrixRankMATLAB(Q,T);
    if(rk==n+2*m) % the system Jacobian is nonsingular
        disp("the system Jacobian is estimated to be nonsingular")
        break
    end
    Ind=double(dual);
    I=zeros(1,m);
    for i=1:numel(Ind)
        if(Ind(i)<=m)
            I(Ind(i))=I(Ind(i))+1;
        elseif(Ind(i)<=2*m)
            I(Ind(i)-m)=I(Ind(i)-m)+1;
        else
            continue
        end
    end
    % row compression and column compression
    [U,S]=compression(Jb(:,I>=2),p);
    [V,T]=compression(Jc(I<2,:).',q);
    V=V.';
    Sd=zeros(n,n);
    Td=zeros(n,n);
    for i =1:n
        Sd(i,i)=S(i);
        Td(i,i)=T(i);
    end
    if rowtrans==true
        [W,~]=compression((Sd*U*J0*V*Td),p);
        U=W*U;
    end
    if rowtrans==false
        [W,~]=compression((Sd*U*J0*V*Td).',q);
        V=V*W;
    end
    % update the matrices
    Upoly=sym(eye(n,n));
    cur_Vpoly=sym(eye(n,n));
    for i=1:n
        for j=1:n
            if(U(i,j)~=0 && p(j)>=p(i))
                Upoly(i,j)=U(i,j)*s^(p(j)-p(i));
            end
            if(V(i,j)~=0 && q(j)>=q(i))
                cur_Vpoly(i,j)=V(i,j)*s^(q(j)-q(i));
            end
        end
    end
    eqs=applyPolynomialMatrixTest.applyPolynomialMatrix(eqs, vars, Upoly);
    % change of variables
    neqs=zeros(n,1,'sym');
    syms x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) x7(t) x8(t) x9(t) x10(t) x11(t) x12(t) x13(t) x14(t) x15(t)
    variables=zeros(n,1,'sym');
    for j=1:n
        variables(j)=str2sym("x"+j+"(t)");
    end
    for i=1:n
        for j=1:n
            if(V(i,j)~=0 && q(j)>=q(i))
                neqs(i)=neqs(i)+V(i,j)*diff(variables(j),q(j)-q(i));
            end
        end
    end
    eqs=subs(eqs,variables,neqs);
end
elapsed=toc

[J,p,q]=systemJacobian(eqs,vars);

J_rand=subs(J,variables,rand(n,1));
J_rand=subs(J_rand,symvar(J_rand),rand(1,numel(symvar(J_rand))));
J_rand=eval(J_rand);
str="after the transformation, the system Jacobian is of size "+num2str(n)+" and of rank "+rank(J);
disp(str);
if(rank(J_rand)~=n)
    error('the transformation failed');
end
