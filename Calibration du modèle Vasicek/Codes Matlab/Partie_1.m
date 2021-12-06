

Y1_plot()
Y2_plot()



    function[] = Y2_plot()
    lim=(0.0075/0.25)-0.5*(0.02/0.25)^2;
        for i=1:100000
            T(i)=i*0.01;
            Y1(i)=Yield(T(i),0.01);
            Y2(i)=Yield(T(i),0.027);
            Y3(i)=Yield(T(i),0.05);
        end
        figure;
        hold on;
        plot(T,Y1);
        plot(T,Y2);
        plot(T,Y3);
        plot(1000,lim,'r*');
        title('Yield en fonction du temps T et selon les valeurs prises par r_0')
        xlabel('T')
        ylabel('Yield')
        legend('r_1=0.01','r_2=0.027','r_3=0.05','limit(Y(r),T,+inf)=0.0268')
        
    end

    function[] = Y1_plot()
        for i=1:3000
            T(i)=i*0.01;
            Y1(i)=Yield(T(i),0.01);
            Y2(i)=Yield(T(i),0.027);
            Y3(i)=Yield(T(i),0.05);
        end
        figure;
        hold on;
        plot(T,Y1);
        plot(T,Y2);
        plot(T,Y3);
        plot(0,0.01,'b*');
        plot(0,0.027,'r*');
        plot(0,0.05,'y*');
        title('Yield en fonction du temps T et selon les valeurs prises par r_0')
        xlabel('T')
        ylabel('Yield')
        legend('r_1=0.01','r_2=0.027','r_3=0.05','limit(Y(r_1),T,0)=r_1','limit(Y(r_2),T,0)=r_2','limit(Y(r_3),T,0)=r_3')
        
    end


    function[B] = B_function(T)
        B = (1-exp(-0.25*T))/0.25;
    end

    function[A] = A_function(T)
        A = (B_function(T)-T)*(0.0075*0.25-0.02^2/2)/0.25^2 - 0.02^2*B_function(T)^2/4*0.25;
    end


    function[Y] = Yield(T,r)
        Y = -(A_function(T)-r*B_function(T))/T;
    end
    

