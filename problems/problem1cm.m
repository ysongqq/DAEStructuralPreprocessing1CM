classdef problem1cm

% PROBLEM1CM    DAE systems that are not contained in the `problem` class of
%               DAEPreprocessingToolbox.

    methods (Static)

        function [eqs, vars] = roboticArmGeneric(K)

        % Robotic arm DAE of size 3*K+2 with symbolic physical parameters.

            syms x0(t)
            syms x(t) [1 K]
            syms phi(t) [1 K]
            syms tau0(t)
            syms tau(t) [1 K]
            syms alp
            syms bet [1 K]
            syms gam [1 K]
            syms Kj [1 K]
            syms Nj [1 K]
            syms J [1 K]
            syms lenth0
            syms lenth [1 K]

            x=x(t);
            phi=phi(t);
            tau=tau(t);
            p0 = symfun(lenth0*cos(1-exp(t)), t);
            p=zeros(K,1,'sym');
            for i=1:K
                p(i) = symfun(lenth0*cos(1-exp(t)) + lenth(i)*cos(1-i*t), t);
            end
            f0=-tau0;
            f=zeros(2*K,1,'sym');
            for i=1:K
                f0=f0-( 2*diff(x0(t))*diff(x(i)) - diff(x(i))^2 )*gam(i)*sin(x(i));
                f(i)=gam(i)*sin(x(i))*diff(x0(t))^2+Kj(i)*(x(i)-phi(i)/Nj(i));
                f(i+K)=(Kj(i)/Nj(i))*(-x(i)+phi(i)/Nj(i))-tau(i);
            end
            Q=alp;
            for i=1:K
                Q=Q+2*gam(i)*cos(x(i))-(bet(i)+gam(i)*cos(x(i)))^2/bet(i);
            end
            Q=1/Q;
            a=zeros(K,1,'sym');
            for i=1:K
                a(i)=-(bet(i)+gam(i)*cos(x(i)))/bet(i);
            end
            eqs0=diff(x0(t),2)+Q*f0;
            for i=1:K
                eqs0=eqs0+Q*a(i)*f(i)-Q*f(i+K);
            end
            eqs1=zeros(K,1,'sym');
            for i=1:K
                eqs1(i)=diff(x(i),2)+Q*a(i)*f0;
                for j=1:K
                    eqs1(i)=eqs1(i)+a(i)*a(j)*Q*f(j);
                    eqs1(i)=eqs1(i)-a(i)*Q*f(j+K);
                    if i==j
                        eqs1(i)=eqs1(i)+f(j)/bet(j);
                    end
                end
            end

            eqs2=zeros(K,1,'sym');
            for i=1:K
                eqs2(i)=diff(phi(i),2)-Q*f0;
                for j=1:K
                    eqs2(i)=eqs2(i)-a(i)*Q*f(j);
                    eqs2(i)=eqs2(i)+Q*f(j+K);
                    if i==j
                        eqs2(i)=eqs2(i)+f(j+K)/J(j);
                    end
                end
            end

            eqs3=cos(x0)-p0;
            eqs4=zeros(K,1,'sym');
            for i=1:K
                eqs4(i)=cos(x0)+cos(x0+x(i))-p(i);
            end
            eqs = [eqs0; eqs1; eqs2; eqs3; eqs4];
            vars = [x0 x phi tau0 tau];
        end

    end
end
