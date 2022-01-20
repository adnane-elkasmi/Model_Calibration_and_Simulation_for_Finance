
%%Partie 1 - L equation de Dupire

function[]=Part_1()
    %% Valeurs initiales définies
    
    %Valeurs utilisées
    Kmax = 20;
    S0 = 10;
    r = 0.1;
    Tmax = 0.5;
    N = 199;
    M = 49;
    
    %Discrétisation de la valeur du strike K et de la maturité T
    K = linspace(0,Kmax,N+2);
    deltaK = Kmax/(N+1);
    t = linspace(0,Tmax,M+2);
    deltat = Tmax/(M+1);

    %Valeurs des paramètres des fonctions
    sigma_VolatiliteFixe = 0.3;
    h = 0.01;
    beta1 = 1;
    beta2 = 1;
    beta = [beta1,beta2];

    
    %% Implémentation des fonctions
    
    %-------------------------------------------------------%
    %-------------Algorithme de Crank-Nicolson--------------%
    %-------------------------------------------------------%
    
    %sigma(i) = 0.3;  % à choisir pour le cas : sigma fixe
    %sigma(i) = beta1/(K(i))^(beta2);     % à choisir pour le cas : sigma variable
    
    function[V] = Prix_Dupire(h,beta1,beta2,sigmaFixeBool)
        
        %Conditions initiales
        for i=1:N+2
            V(1,i) = max(S0-K(i),0);
        end

        %Sigma (variable ou fixe)
        if sigmaFixeBool == 1
            %Si on choisit la volatilité fixe
            for i = 1:N+2
                sigma(i) = sigma_VolatiliteFixe + h; 
            end
        else
            %Si on choisit la volatilité variable
            for i = 1:N+2
                sigma(i) = (beta1/((K(i))^(beta2))) + h; 
            end
        end
        


        %Conditions aux limites
        for n = 1:M+2
            V(n,1) = S0;
            V(n,N+2) = 0;
        end

        %Calcul des Ai, Bi et Di
        for i = 2:N+1
            B(i) = -(deltat/4)*(K(i)*r/deltaK + ((sigma(i)*K(i))/deltaK)^2);
            A(i) = (deltat/4)*(K(i)*r/deltaK - ((sigma(i)*K(i))/deltaK)^2);
            D(i) = 1 + (deltat/2)*((sigma(i)*K(i)/deltaK)^2) ;
        end



        %Calcul des Ci
        for n = 1:M+1
            
            for i = 2:N+1
                if i == 2
                    C(n,i)=-(deltat/4)*(K(i)*r/deltaK - (sigma(i)^2 * K(i)^2)/deltaK^2)*V(n,i+1) + (1-deltat/2*sigma(i)^2*K(i)^2/deltaK^2)*V(n,i) + deltat/4*(r*K(i)/deltaK + sigma(i)^2 * K(i)^2 /deltaK^2)*V(n,i-1) - S0*B(2);
                else
                    C(n,i)=-(deltat/4)*(K(i)*r/deltaK - (sigma(i)^2 * K(i)^2)/deltaK^2)*V(n,i+1) + (1-deltat/2*sigma(i)^2*K(i)^2/deltaK^2)*V(n,i) + deltat/4*(r*K(i)/deltaK + sigma(i)^2 * K(i)^2 /deltaK^2)*V(n,i-1);
                end
            end
            
            C_star(n,2) = C(n,2);
            D_star(2) = D(2);

            for i=3:N+1
                D_star(i)=D(i)-(B(i)*A(i-1))/D_star(i-1);
                C_star(n,i)=C(n,i)-(B(i)*C_star(n,i-1))/D_star(i-1);
            end

            %Calcul de V
            V(n+1,N+1)=C_star(n,N+1)/D_star(N+1);

            for i=N:-1:2
                V(n+1,i)=(C_star(n,i)-A(i)*V(n+1,i+1))/D_star(i);
            end
            
        end

        
    end


    %----------------------------------------------------------%
    %----Fonction permettant de calculer la Vega de l'option---%
    %----------------------------------------------------------%
    function[vega]=Vega_Dupire(h,beta1,beta2,sigmaFixeBool)
        vega=(Prix_Dupire(h,beta1,beta2,sigmaFixeBool)-Prix_Dupire(0,beta1,beta2,sigmaFixeBool))/h;
    end



    %% Tests réalisés

    %Résoudre l'équation de Dupire (OK) et trouver la matrice V(n,i) pour sigma
    %fixe 
    V_sigmaFixe = Prix_Dupire(0,beta1,beta2,1);
    disp(V_sigmaFixe);
    
    %Résoudre l'équation de Dupire (OK) et trouver la matrice V(n,i) pour sigma
    %variable 
    V_sigmaVariable = Prix_Dupire(0,beta1,beta2,0);
    disp(V_sigmaVariable);

    
    
    %Visualisation de la fonction V pour les maturités T, T/2 et T=0 en 2D et
    %3D pour sigma fixe 
        figure;
        plot(K,V_sigmaFixe(1,:),K,V_sigmaFixe(floor((M+1)/2),:),K,V_sigmaFixe(M+1,:));
        legend('t=0','t=T/2','t=T','location','northeast');
        title('Graph de V en fonction de K et selon T - volatilité fixe')
        xlabel('K');
        ylabel('V');

        %en 3D
        figure;
        mesh(K,t,V_sigmaFixe);
        title('Surface de V selon T et K - volatilité fixe');
        xlabel('K');
        ylabel('T');
        zlabel('V');

    %Visualisation de la fonction V pour les maturités T, T/2 et T=0 en 2D et
    %3D pour sigma variable
        figure;
        plot(K,V_sigmaVariable(1,:),K,V_sigmaVariable(floor((M+1)/2),:),K,V_sigmaVariable(M+1,:));
        legend('t=0','t=T/2','t=T','location','northeast');
        title('Graph de V en fonction de K et selon T - volatilité variable')
        xlabel('K');
        ylabel('V');

        %en 3D
        figure;
        mesh(K,t,V_sigmaVariable);
        title('Surface de V selon T et K - volatilité variable');
        xlabel('K');
        ylabel('T');
        zlabel('V');

        
        
    %Determination de la Vega de l'option pour un sigma Fixe
    Vega_sigmaFixe = Vega_Dupire(h,beta1,beta2,1);
    
    %Determination de la Vega de l'option pour un sigma Variable
    Vega_sigmaVariable = Vega_Dupire(h,beta1,beta2,0);
    
    
    
    %Visualisation de la fonction Vega en 2D et 3D pour sigma fixe
        figure;
        plot(K,Vega_sigmaFixe(1,:),K,Vega_sigmaFixe(floor((M+1)/2),:),K,Vega_sigmaFixe(floor(M+2),:));
        legend({'t=0','t=T/2','t=T','location','northeast'});
        xlabel('K');
        ylabel('Vega');
        title('Graph de Vega en fonction de K et selon T - volatilité fixe');

        %en 3D
        figure;
        mesh(K,t,Vega_sigmaFixe);
        title('Surface de Vega en fonction de T et K - volatilité fixe');
        xlabel('K');
        ylabel('T');
        zlabel('Vega');
        
    %Visualisation de la fonction Vega en 2D et 3D pour sigma variable
        figure;
        plot(K,Vega_sigmaVariable(1,:),K,Vega_sigmaVariable(floor((M+1)/2),:),K,Vega_sigmaVariable(floor(M+2),:));
        legend({'t=0','t=T/2','t=T','location','northeast'});
        xlabel('K');
        ylabel('Vega');
        title('Graph de Vega en fonction de K et selon T - volatilité variable');

        %en 3D
        figure;
        mesh(K,t,Vega_sigmaVariable);
        title('Surface de Vega en fonction de T et K - volatilité variable');
        xlabel('K');
        ylabel('T');
        zlabel('Vega');


end
