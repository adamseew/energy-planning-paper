
%%
% load the data

% change the following 2
load('new_data_to_sort/3casadi_4.mat'); % load the data
filename = '3casadi_4';

n_coefs = 2*r+1;
titles_coefs = {'\alpha_0','\alpha_1','\beta_1','\alpha_2','\beta_2',...
    '\alpha_3','\beta_3','\alpha_4','\beta_4'}; % max r=4

x0=10;
y0=10;
width=1500;
height=1000;


%% path
h1 = figure('DefaultTextFontSize',24,'Color','w');
set(gcf,'position',[x0,y0,width,width]);
   
path_fn = strcat('path_',filename,'.gif');

for j=1:1/delta_T:k % path
    plot(log_p(1,1:j),log_p(2,1:j),'Color','k','LineWidth',2)
    grid on
    ax = gca;
    ax.FontName = 'Helvetica';
    ax.FontSize = 24;
    xlim([-150 200]);
    xlabel('x (m)','FontName','Helvetica');
    ylim([-100 250]);
    ylabel('y (m)','FontName','Helvetica');
    
    frame = getframe(h1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,32);
    
    if j == 1 
        imwrite(imind,cm,path_fn,'gif','Loopcount',1); 
    else 
        imwrite(imind,cm,path_fn,'gif','WriteMode','append'); 
    end 
end


%% coefs

h2 = figure('DefaultTextFontSize',24,'Color','w');
set(gcf,'position',[x0,y0,width/2,n_coefs*height/4]);
   
coefs_fn = strcat('coefs_',filename,'.gif');

for j=1:1/delta_T:k % path
    
    t = tiledlayout(n_coefs,1);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    for jj=1:n_coefs
        coef_data = log_q(jj,1:j);
        nexttile,plot(linspace(0,j*delta_T,length(coef_data)),coef_data,...
            'Color','k','LineWidth',1.2)
        grid on
        xlim([0 k*delta_T]);
        ax = gca;
        ax.FontName = 'Helvetica';
        ax.FontSize = 24;
        title(string(titles_coefs(jj)),'Interpreter','tex');
    end
    xlabel(t,'Time (sec)','FontName','Helvetica','FontSize',24)
    ylabel(t,'Value','FontName','Helvetica','FontSize',24)
    
    frame = getframe(h2);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,32);
    
    if j == 1 
        imwrite(imind,cm,coefs_fn,'gif','Loopcount',1); 
    else 
        imwrite(imind,cm,coefs_fn,'gif','WriteMode','append'); 
    end 
end


%% ener

h3 = figure('DefaultTextFontSize',24,'Color','w');
set(gcf,'position',[x0,y0,width,height]);
   
ener_fn = strcat('ener_',filename,'.gif');

data_mpc_q = nan;

for j=1:1/delta_T:k % path
    
    if (j-1)*delta_T >= 85
        data_mpc_q = log_q_chain(log_q_chain(:,1) == (j-1)*delta_T,:);
        data_mpc_q = data_mpc_q(2:end);
        data_mpc_q = [linspace((j-1)*delta_T,(j-1+N*(1/delta_T))*...
            delta_T,N).' data_mpc_q.'];
    end
    
    if ~isnan(data_mpc_q)
        plot(data_mpc_q(:,1),data_mpc_q(:,2),'Color','r','LineWidth',1.2)
    end
    
    ener_data = log_y(1:j);
    hold on
    plot(log_b(:,1),log_b(:,3),'Color','g','LineWidth',1.2)
    plot(linspace(0,j*delta_T,length(ener_data)),ener_data,...
        'Color','k','LineWidth',1.2)
    hold off
    
    grid on
    xlim([0 k*delta_T]);
    ylim([25 50]);
    ax = gca;
    ax.FontName = 'Helvetica';
    ax.FontSize = 24;
    xlabel('Time (sec)','FontName','Helvetica','FontSize',24)
    ylabel('Power (W)','FontName','Helvetica','FontSize',24)
    
    frame = getframe(h3);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,32);
    
    if j == 1 
        imwrite(imind,cm,ener_fn,'gif','Loopcount',1); 
    else 
        imwrite(imind,cm,ener_fn,'gif','WriteMode','append'); 
    end 
end

%% ctls

h4 = figure('DefaultTextFontSize',24,'Color','w');
set(gcf,'position',[x0,y0,width,height/2]);
   
ener_fn = strcat('ctls_',filename,'.gif');

data_mpc_c = nan;

for j=1:1/delta_T:k % path
    
    if (j-1)*delta_T >= 85
        data_mpc_c = log_c_chain(log_c_chain(:,1) == (j-1)*delta_T,:);
        data_mpc_c = data_mpc_c(2:end);
        data_mpc_c = imresize(data_mpc_c,[1 (1/delta_T)*(N-1)]);
        data_mpc_c = [linspace((j-1)*delta_T,(j-1+(N-1)*(1/delta_T))*...
            delta_T,length(data_mpc_c)).' data_mpc_c.'];
    end
    
    
    t = tiledlayout(1,2);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    c1_data = log_c1(1:j);
    c2_data = log_c2(1:j);
    nexttile,plot(linspace(0,j*delta_T,length(c1_data)),c1_data,...
        'Color','k','LineWidth',1.2)
    grid on
    ylim([-1000 0]);
    xlim([0 k*delta_T]);
    ax = gca;
    ax.FontName = 'Helvetica';
    ax.FontSize = 24;
    title('c_{i,1}','Interpreter','tex');
    nexttile
    if ~isnan(data_mpc_c)
        plot(data_mpc_c(:,1),data_mpc_c(:,2),'Color','r','LineWidth',1.2)
    end
    hold on
    plot(linspace(0,j*delta_T,length(c2_data)),c2_data,...
        'Color','k','LineWidth',1.2)
    hold off
    grid on
    ylim([2 10]);
    xlim([0 k*delta_T]);
    ax = gca;
    ax.FontName = 'Helvetica';
    ax.FontSize = 24;
    title('c_{i,2}','Interpreter','tex');
    
    
    xlabel(t,'Time (sec)','FontName','Helvetica','FontSize',24)

    
    frame = getframe(h4);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,32);
    
    if j == 1 
        imwrite(imind,cm,ener_fn,'gif','Loopcount',1); 
    else 
        imwrite(imind,cm,ener_fn,'gif','WriteMode','append'); 
    end 
end



