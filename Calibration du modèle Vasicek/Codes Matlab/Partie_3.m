%%TP3 - Vasicek Model
%%Part 3 - Levenberg Marquart

function[]=Partie_3()
    %% Valeurs initiales

    epsilon = 10^-9;
    r0 = 0.04;
    %r0 = 0.025;
    lambda = 0.01;
    

    %% Fonctions préalablement définies

function[B] = B_function(T,gamma)
    B = (1-exp(-gamma*T))/gamma;
end

function[A] = A_function(T,eta,sigmacarre,gamma)
    B = B_function(T,gamma);
    A = (B-T)*(eta*gamma-sigmacarre/2)/gamma^2-sigmacarre*B^2/(4*gamma);
end


function[Y] = Yield_function(T, eta, sigmacarre, gamma, r0)
    Y = -(A_function(T,eta,sigmacarre,gamma)-r0*B_function(T,gamma))/T;
end

function[dB] = dB_function(t,T,gamma)
    Tau = T - t;
    B = B_function(Tau,gamma);
    dB = (Tau * exp(-gamma*Tau)-B)/gamma;
end

 function[dA] = dA_function(t,T,Beta)
     Tau = T - t;
     eta = Beta(1);
     sigmacarre = Beta(2);
     gamma = Beta(3);
     B = B_function(Tau, gamma);
     dB = dB_function(t,T,gamma);
     dA = (1/gamma^2)*(eta*(dB*gamma-B)+Tau*eta-(sigmacarre)/2*(dB-2*B/gamma)-Tau*(sigmacarre)/gamma - (sigmacarre)*B/4*(2*gamma*dB-B));
 end

    %% Définition de la fonction pour la calibration de la courbe des taux d'intérêts à t

        
    function[]=ExecutionCalibration(t,r0,lambda,epsilon,T,YieldMarket)
        
        %Initialisation
        J = zeros(10,3);
        Rth = zeros(10,1);
        Res=zeros(10,1);

        eta = 1;
        sigma = 1 ; 
        gamma = 1;
        Beta = [eta,sigma^2,gamma];
    
        d = [1,1,1];    
        
        k = 0;

        %Calibration
        while norm(d) > epsilon 
            
            eta = Beta(1);
            sigmacarre = Beta(2);
            gamma = Beta(3);
            
            for p=1:10
                
               Tau = T(p)-t;
               A = A_function(Tau, eta, sigmacarre, gamma);
               B = B_function(Tau,gamma);
               dA = dA_function(t,T(p),Beta);
               dB = dB_function(t, T(p),gamma);
               
               J(p,1)=(B-Tau)/(Tau*eta);
               J(p,2)=(-1/(Tau*gamma))*((B-Tau)/(2*gamma)+(B^2)/4);
               J(p,3) = (1/Tau)*(dA-r0*dB);
               Res(p) = YieldMarket(p)+(A-r0*B)/Tau;
               
            end

            TransposeJ = J.';
            M = TransposeJ*J + lambda*eye(3);
            d = -inv(M)*TransposeJ*Res;
            Beta = Beta+d.';
            
            k = k+1;
            
        end
        
        %Construction du vecteur Yiel
        for p=1:10
            
            eta = Beta(1);
            sigmacarre = Beta(2);
            gamma = Beta(3);
            Tau = T(p)-t;
            Rth(p) = Yield_function(Tau, eta, sigmacarre, gamma, r0);
            
        end

        %%Tracé des figures
        grid on;
        hold on;
        plot (T, YieldMarket,'*'); %Valeurs du marché
        plot(T,Rth,'r'); %Tracé de la courbe des taux
        hold off;
        title(['Courbe des taux dintérêts calibrée pour p = 10 valeurs et à t = ', num2str(t)]);
        legend({'Valeurs du marché : taux zéro-coupon R','Courbe des taux théoriques Rth'},'Location','southeast','Orientation','vertical')
        xlabel('Maturité T(p)') 
        ylabel('Taux YieldMarket(p)')

        %%Valeurs des Beta et de k
        disp(['Affichage des différentes valeurs du vecteur Beta pour t = ',num2str(t)]) ;
        disp(['La valeur de eta (après calibration) est ', num2str(Beta(1))]) ;
        disp(['La valeur de sigmacarre (après calibration) est ', num2str(Beta(2))]) ;
        disp(['La valeur de gamma (après calibration) est ', num2str(Beta(3))]) ;
        disp(['La valeur de k est ', num2str(k)]);
        disp('--------------------------------------------');
        
    end

    %% Première execution avec t = 0
    
    YieldMarket=[0.035,0.041,0.0439,0.046,0.0484,0.0494,0.0507,0.0514,0.052,0.0523];
    T=[3,6,9,12,15,18,21,24,27,30];
    
    ExecutionCalibration(0,r0,lambda,epsilon,T,YieldMarket)


    %% Seconde execution avec t = 1
    
    YieldMarket=[0.056,0.064,0.074,0.081,0.082,0.09,0.087,0.092,0.0895,0.091];
    T=[3,6,9,12,15,18,21,24,27,30];
    
    figure;
    
    ExecutionCalibration(1,r0,lambda,epsilon,T,YieldMarket)

end

