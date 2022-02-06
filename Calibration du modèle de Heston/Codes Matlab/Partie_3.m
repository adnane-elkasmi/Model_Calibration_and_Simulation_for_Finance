%TP5 - Calibration du modèle de Heston

%Partie 3 - Smile de Volatilité

function[] = Partie_3()
    %% Valeurs initiales
    
    S0=10;
    T=0.5;
    v0=0.03;
    N=100;
    deltat=T/N;
    r=0.1;
    k=0.3;
    rho=0.7;
    theta=0.3;
    eta=0.4;
    eps=10^-3;

    %Le vecteur strike prend les valeurs de 5 à 20 avec un pas de 0.5)
    K=linspace(5,20,31);
    
    %Nombre de simulations Monte-Carlo
    Nmc=10000;
    

    %% Implémentations des fonctions principales

    %Calcul de l'actif S et de la volatilité v
    function[S,v,S_sym,v_sym]=Calcul_ActifVola()

        S(1)=S0;
        v(1)=v0;
        S_sym(1)=S0;
        v_sym(1)=v0;

        for i=1:N
            N1=randn;
            N2=randn;
            v(i+1)=v(i)+k*(theta-v(i))*deltat+eta*sqrt(abs(v(i)))*sqrt(deltat)*N1+(eta^2/4)*(deltat*N1^2-deltat);
            S(i+1)=S(i)*exp((r-v(i)/2)*deltat+sqrt(abs(v(i)))*(rho*sqrt(deltat)*N1+sqrt(1-rho^2)*sqrt(deltat)*N2));
            v_sym(i+1)=v_sym(i)+k*(theta-v_sym(i))*deltat+eta*sqrt(abs(v_sym(i)))*sqrt(deltat)*(-N1)+(eta^2/4)*(deltat*(-N1)^2-deltat);
            S_sym(i+1)=S_sym(i)*exp((r-v_sym(i)/2)*deltat+sqrt(abs(v_sym(i)))*(rho*sqrt(deltat)*(-N1)+sqrt(1-rho^2)*sqrt(deltat)*(-N2)));
        end
        
    end


    %Calcul du pay off du call européen
    function[f] = Payoff_Call_Europeen(S,K)
        f = max(S-K, 0);
    end


    %Calcul du prix du call selon l'estimateur 2
    function[Prix] = Prix_Call_Europeen_S0_fixe_Estimateur2(K)
        Somme = 0;
        
        for n = 1 : Nmc
            [S,~,S_sym,~]=Calcul_ActifVola();
            Somme = Somme + Payoff_Call_Europeen(S(N+1),K) + Payoff_Call_Europeen(S_sym(N+1),K);
        end
        
        Prix = Somme/(2*Nmc);
    end
    
    function[V_smile]=Simul_V(S0,K,r,T)
        for j = 1 : 31
            Prix=Prix_Call_Europeen_S0_fixe_Estimateur2(K(j));
            
            while max(S0-K(j)*exp(-r*T),0)>Prix || Prix>S0
                Prix=Prix_Call_Europeen_S0_fixe_Estimateur2(K(j));
            end
            
            V_smile(j)=Prix;
        end
    end


    %Simulation de V
    V = Simul_V(S0,K,r,T);
    
    
    %Fonction de répartition de la loi Normale
    function[g]=Nx(x)
        g=0.5*(1+erf(x/sqrt(2)));
    end

    %Fonctions d1,d2 VBS et dVBS
    function[d1]=d1_function(S0,K,r,sigma,T)
        d1=(log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
    end


    function[d2]=d2_function(S0,K,r,sigma,T)
        d2=(log(S0/K)+(r-sigma^2/2)*T)/(sigma*sqrt(T));
    end

    %Fonction du Call théorique
    function[f]=Vbs(S0,K,r,sigma,T)
        f=S0*Nx(d1_function(S0,K,r,sigma,T))-K*exp(-r*T)*Nx(d2_function(S0,K,r,sigma,T));
    end


    function[F]=dVbs(S0,K,r,sigma,T)
        F=S0*sqrt(T/(2*pi))*exp(-d1_function(S0,K,r,sigma,T)^2/2);
    end
    
    
    %Algorithme de Newton
    function [sigmaNewton]=Algo_Newton(k)
        
        sigmaNewton=sqrt(2*abs((log(S0/K(k))+r*T)/T));
        
        while Vbs(S0,K(k),r,sigmaNewton,T)-V(k)>eps
            sigmaNewton=sigmaNewton-(Vbs(S0,K(k),r,sigmaNewton,T)-V(k))/dVbs(S0,K(k),r,sigmaNewton,T);
        end
        
    end

   

    %% Tracé du smile de volatilité
    
    function[] = Smile_Volatilite()
        

        for l=1:31
            sigmaImplicite(l)=Algo_Newton(l);
            PrixBS(l)=Vbs(S0,K(l),r,sigmaImplicite(l),T);
        end

        
        %Tracé des figures
        figure;
        plot(K,sigmaImplicite);
        title('Volatilité implicite en fonction de K')
        xlabel('K');
        ylabel('Sigma_{imp}');
        
        figure;
        plot(K,V,'o',K,PrixBS);
        legend('V_{Heston}','V_{BS}');
        title('Prix Heston et Black-Scholes en fonction de K')
        xlabel('K');
        ylabel('Prix');
        
    end


    %% Appel des fonctions / Tests
    
    Smile_Volatilite();
    
    
end
