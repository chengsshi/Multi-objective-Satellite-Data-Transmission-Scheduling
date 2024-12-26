function pareto_plot(trace_obj, trace_con, maxgen)

    figure('Position', [100 100 640 480]);
    
    for i = 1:maxgen
  
        hold on;
      
        plot(trace_obj(maxgen-i+1), trace_con(i), 'bo');
     
        grid on;
   
        xlabel('传输时间');
        ylabel('卫星和地面站的优先级乘积');
    
        title('Pareto图');
        
 
        legend('Pareto最优解');
        
   
        hold off;
        
 
        drawnow;
    end
end


