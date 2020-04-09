clear
clc
close all
% basic parameter




f = figure(1);

for seq = 1:4
    
    j = 1;
    lambda = linspace(750,1500,2000) * 1E-9;
    theta = pi / 3;
    n_core = 1.444;
    n_clad = 1;
    p_core = n_core^2;
    p_clad = n_clad^2;
    metal_t = 75E-9;
    n_glass = 1.22;
    p_glass = n_glass^2;
%     L = 50E-6;
%     D = 2E-6;
    L_div_D = 25;
    

    switch seq
        case 1
            m = linspace(30,105,6); % thinckness
            mark_thinckness = {[strcat('d_m=', num2str(m(1))), 'nm'], [strcat('d_m=', num2str(m(2))), 'nm'],...
                            [strcat('d_m=', num2str(m(3))), 'nm'], [strcat('d_m=', num2str(m(4))), 'nm'],...
                            [strcat('d_m=', num2str(m(5))), 'nm'], [strcat('d_m=', num2str(m(6))), 'nm']};
        case 2
            m = linspace(20,70,6);   % L / D
            mark_L_div_D = {strcat('L/D=', num2str(m(1))), strcat('L/D=', num2str(m(2))),...
                            strcat('L/D=', num2str(m(3))), strcat('L/D=', num2str(m(4))),...
                            strcat('L/D=', num2str(m(5))), strcat('L/D=', num2str(m(6)))};
        case 3
            m = linspace(1.44,1.4436,6);  % n_core
            n_core_max = max(m);
            n_core_min = min(m);
            n_core_range = m;
            mark_n_core = {strcat('n_{core}=', num2str(m(1), '%.4f')),strcat('n_{core}=', num2str(m(2), '%.4f')),...
                            strcat('n_{core}=', num2str(m(3), '%.4f')),strcat('n_{core}=', num2str(m(4), '%.4f')),...
                            strcat('n_{core}=', num2str(m(5), '%.4f')),strcat('n_{core}=', num2str(m(6), '%.4f'))};
%                             strcat('n_{core}=', num2str(m(7), '%.4f')),strcat('n_{core}=', num2str(m(8), '%.4f'))};
        case 4
            m = linspace(1.22,1.2232,6); % n_glass
            n_glass_max = max(m);
            n_glass_min = min(m);
            n_glass_range = m;
            mark_n_glass = {strcat('n_d=', num2str(m(1), '%.4f')),strcat('n_d=', num2str(m(2), '%.4f')),...
                            strcat('n_d=', num2str(m(3), '%.4f')),strcat('n_d=', num2str(m(4), '%.4f')),...
                            strcat('n_d=', num2str(m(5), '%.4f')),strcat('n_d=', num2str(m(6), '%.4f'))};
%                             strcat('n_d=', num2str(m(7), '%.4f')),strcat('n_d=', num2str(m(8), '%.4f'))};
    end
            
    
    for d = m
        
        switch seq
            case 1
                metal_t = d * 1E-9;
            case 2
                L_div_D = d;
            case 3
                n_core = d;
                p_core = n_core^2;
            case 4
                n_glass = d;
                p_glass = n_glass^2;
        end

%         metal_t = d * 1E-9;

        q1 = sqrt(p_core - n_core^2 * (sin(theta))^2) / p_core;
        % q1 = sqrt(1 / p_core) * cos(theta);

    %     n_glass = d;
    %     p_glass = n_glass^2;

        p_metal = Au_epsilon(lambda);

        % refractive
        n_metal = sqrt(p_metal);

        % layer 2: metal

        theta2 = asin(sin(theta) ./ (n_core ./ n_metal));
        q2 = sqrt(p_metal - n_core^2 .* (sin(theta))^2) ./ p_metal;
    %     q2 = sqrt(1 / p_metal) * cos(theta2);
        beta2 = 2 * pi .* metal_t ./ lambda .* sqrt(p_metal - n_core^2 * sin(theta)^2);

        % layer 3: glass
        theta3 = asin(sin(theta2) ./ (n_metal ./ n_glass));
        q3 = sqrt(p_glass - n_core^2 * (sin(theta))^2) / p_glass;
    %     q3 = sqrt(1 / p_glass) * cos(theta3);

        % transfer matrix

        m11 = cos(beta2);
        m12 = -1i .* sin(beta2) ./ q2;
        m21 = -1i .* q2 .* sin(beta2);
        m22 = m11;

        % power

        p_trans = n_core^2 * sin(theta) * cos(theta) / (1 - n_core^2 * cos(theta)^2)^2;

        % reflection ratio

        r_p = ((m11 + m12 .* q3) .* q1 - (m21 + m22 .* q3)) ./ ((m11 + m12 .* q3) .* q1 + (m21 + m22 .* q3));
        R_p = abs(r_p).^2;
        Nef = floor(L_div_D / tan(theta) / 2);
        if(seq == 2)
            p = real(10*log(R_p.^Nef));
        else
            p = R_p.^Nef;
        end
        switch seq
            case 3
                [min_y_1, min_x_1(j)] = min(p);
            case 4
                [min_y_2, min_x_2(j)] = min(p);
        end
        j = j + 1;
        switch seq
            case 1
                 fig_1 = subplot(2, 2, seq);
            case 2
                 fig_2 = subplot(2, 2, seq);
            case 3
                 fig_3 = subplot(2, 2, seq);
            case 4
                 fig_4 = subplot(2, 2, seq);
        end
%         fig_1 = subplot(2, 2, seq);
        plot(lambda*1E9, p, 'LineWidth', 1);
        hold on;
        xlabel('Wavelength(nm)','FontSize', 12, 'LineWidth', 1);
        if(seq == 2)
            ylabel('Normalization Transmisson(dB)','FontSize', 12, 'LineWidth', 1);
            axis([min(lambda)*1E9 max(lambda)*1E9 -120 0]);
        else
            ylabel('Normalization Transmisson','FontSize', 12, 'LineWidth', 1);
            axis([min(lambda)*1E9 max(lambda)*1E9 0 1]);
        end
%         ylabel('Normalization Transmisson','FontSize', 12, 'LineWidth', 1);
%         axis([min(lambda)*1E9 max(lambda)*1E9 0 1]);
        legend('boxoff');
    end
end

legend(fig_1, mark_thinckness, 'NumColumns', 3, 'location', 'north');

legend(fig_2, mark_L_div_D, 'NumColumns', 2, 'location', 'southeast');

legend(fig_3, mark_n_core, 'NumColumns', 3, 'location', 'north');

legend(fig_4, mark_n_glass, 'NumColumns', 3, 'location', 'north');

dim = [0.05 0.7 0.3 0.3];
h = annotation('textbox',dim,'String','(a)','FitBoxToText','on','FontSize', 16);
h.LineStyle = 'none';

dim = [0.54 0.7 0.3 0.3];
h = annotation('textbox',dim,'String','(b)','FitBoxToText','on','FontSize', 16);
h.LineStyle = 'none';

dim = [0.05 0.2 0.3 0.3];
h = annotation('textbox',dim,'String','(c)','FitBoxToText','on','FontSize', 16);
h.LineStyle = 'none';

dim = [0.54 0.2 0.3 0.3];
h = annotation('textbox',dim,'String','(d)','FitBoxToText','on','FontSize', 16);
h.LineStyle = 'none';


f_2 = figure(2);

f_2_1 = subplot(1,2,1);
x = linspace(min(n_core_min), max(n_core_max), 1000);
val_1 = polyfit(n_core_range, lambda(min_x_1)*1e9, 1);
y = polyval(val_1, x);

plot(n_core_range, lambda(min_x_1)*1e9, 'o', 'LineWidth', 1);
hold on;
plot(x, y, '--', 'LineWidth', 1);
axis([min(n_core_min) max(n_core_max) min(lambda(min_x_1))*1e9 max(lambda(min_x_1))*1e9 + 10]);
xlabel('n_{core}','FontSize', 12, 'LineWidth', 1);
ylabel('Wavelength(nm)','FontSize', 12, 'LineWidth', 1);
legend('Origin data', 'Fit curve');

f_2_2 = subplot(1,2,2);
x = linspace(min(n_glass_min), max(n_glass_max), 1000);
val_2 = polyfit(n_glass_range, lambda(min_x_2)*1e9, 1);
y = polyval(val_2, x);

plot(n_glass_range, lambda(min_x_2)*1e9, 'o', 'LineWidth', 1);
hold on;
plot(x, y, '--', 'LineWidth', 1);
axis([min(n_glass_min) max(n_glass_max) min(lambda(min_x_2))*1e9 max(lambda(min_x_2))*1e9 + 10]);
xlabel('n_d','FontSize', 12, 'LineWidth', 1);
ylabel('Wavelength(nm)','FontSize', 12, 'LineWidth', 1);
legend('Origin data', 'Fit curve');
dim = [0.04 0.7 0.3 0.3];
h = annotation('textbox',dim,'String','(a)','FitBoxToText','on','FontSize', 16);
h.LineStyle = 'none';
dim = [0.25 0.25 0.3 0.3];
h = annotation('textbox',dim,'String',['y=',num2str(val_1(1)),'x+',num2str(val_1(2))],'FitBoxToText','on','FontSize', 10);
h.LineStyle = 'none';

dim = [0.54 0.7 0.3 0.3];
h = annotation('textbox',dim,'String','(b)','FitBoxToText','on','FontSize', 16);
h.LineStyle = 'none';
dim = [0.75 0.14 0.3 0.3];
h = annotation('textbox',dim,'String',['y=',num2str(val_2(1)),'x',num2str(val_2(2))],'FitBoxToText','on','FontSize', 10);
h.LineStyle = 'none';

figWidth = 10; % 设置图片宽度

figHeight = 8;  % 设置图片高度

set(f,'PaperUnits','inches'); % 图片尺寸所用单位

set(f,'PaperPosition',[0 0 figWidth figHeight]);

set(f_2,'PaperUnits','inches'); % 图片尺寸所用单位

set(f_2,'PaperPosition',[0 0 10 4]);

set(f_2_1,'Position',[0.07 0.14 0.4 0.76]);
set(f_2_1,'FontName','Helvetica','FontSize', 12, 'LineWidth', 1);
set(f_2_2,'Position',[0.57 0.14 0.4 0.76]);
set(f_2_2,'FontName','Helvetica','FontSize', 12, 'LineWidth', 1);


set(fig_1,'FontName','Helvetica','FontSize', 12, 'LineWidth', 1);

set(fig_1,'Position',[0.06 0.57 0.42 0.36]);


set(fig_2,'FontName','Helvetica','FontSize', 12, 'LineWidth', 1);

set(fig_2,'Position',[0.56 0.57 0.42 0.36]);


set(fig_3,'FontName','Helvetica','FontSize', 12, 'LineWidth', 1);

set(fig_3,'Position',[0.06 0.07 0.42 0.36]);


set(fig_4,'FontName','Helvetica','FontSize', 12, 'LineWidth', 1);

set(fig_4,'Position',[0.56 0.07 0.42 0.36]);


print(f,'four_para.tif','-r600','-dtiff'); % 设置图片格式、分辨率
print(f_2,'senstivity_para.tif','-r600','-dtiff'); % 设置图片格式、分辨率