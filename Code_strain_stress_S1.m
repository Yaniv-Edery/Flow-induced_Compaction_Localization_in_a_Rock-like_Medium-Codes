%% %%%%%%%%% PREPARING THE DATA %%%%%%%%%%%
clc
clear
close all
datapath='E:\Compaction_localization_data\S1';
load(fullfile(datapath, 'Workspace_strain_stress_S1_updated.mat'));
% load(fullfile(datapath, 'Workspace_full_S1.mat'),'u_original','v_original');
% load(fullfile(datapath, 'Workspace_press_disch_analyzed_S1.mat'),'P_new','Q_new');

% Orientation
u1={};
for i=1:length(u_original)
    u1{i}=flipud(fliplr(-u_original{i}));
end

u2={};
for i=1:length(v_original)
   u2{i}=flipud(fliplr(-v_original{i}));
end

% Dimensions
L=0.029; %length in [m]
H=0.0085; %height in [m]
sz=size(u1{1});

%% %%%%% CALCULATING THE FLUX [cm/s] %%%%%%%%
peak=find(P_new==max(P_new)); %[Pa]
q=Q_new./(H*10^2*0.7); %[cm/s]
figure
hold on
plot(P_new(1:peak).*10^-6,q(1:peak),'Linewidth',1.5)
plot(P_new(peak+1:end).*10^-6,q(peak+1:end),'Linewidth',1.5)
legend('Increase','Decrease','location','northwest')
xlabel('Pressure [MPa]')
ylabel('Flux [cm/s]')
xlim([0 0.7])
ylim([0 1.5])
set(gca,'FontSize',18)
grid on

%calculating the maximum Reynolds number 
ro=830; %Fluid's density [kg/m^3]
d=96*10^-6; %Mean grain size in [m]
mu=0.0117; %Fluid's viscosity[Pa.s]
R_max=(ro*(max(q)*10^-2)*d)/mu

%% %%%%%% SMOOTHING %%%%%%%%%%%

window_size = 6; % Smoothing window size (e.g., 5 for 5x5 neighborhood)
boundary_width = 2; % Number of boundary cells to set as NaN (adjustable)
num_cells = length(u1); 

% Creating empty cell arrays
u1_affine = cell(size(u1)); 
u2_affine = cell(size(u2));


for c = 1:num_cells
    % Extract displacement fields from current cell
    U1 = u1{c};
    U2 = u2{c};
    
    % Get grid size
    [N2,N1] = size(U1);  
    
    % Apply NaN boundary mask
    U1(1:boundary_width, :) = NaN;
    U1(end-boundary_width+1:end, :) = NaN;
    U1(:, 1:boundary_width) = NaN;
    U1(:, end-boundary_width+1:end) = NaN;
    
    U2(1:boundary_width, :) = NaN;
    U2(end-boundary_width+1:end, :) = NaN;
    U2(:, 1:boundary_width) = NaN;
    U2(:, end-boundary_width+1:end) = NaN;
    
    % Creating empty displacement fields
    U1_smooth = NaN(size(U1));
    U2_smooth = NaN(size(U2));
    
    half_w = floor(window_size / 2); % Half window size for indexing
    
    for i = 1:N2
        for j = 1:N1
            % Skip computation if the current point is NaN
            if isnan(U1(i, j)) || isnan(U2(i, j))
                continue;
            end

            % Define local neighborhood bounds making sure it is not out of
            % the domain
            i_min = max(i - half_w, 1);
            i_max = min(i + half_w, N2);
            j_min = max(j - half_w, 1);
            j_max = min(j + half_w, N1);
            
            % Get local coordinates
            [X, Y] = meshgrid(j_min:j_max, i_min:i_max);
            X = X(:);
            Y = Y(:);
            
            % Get local displacements
            U1_local = U1(i_min:i_max, j_min:j_max);
            U2_local = U2(i_min:i_max, j_min:j_max);
            U1_local = U1_local(:);
            U2_local = U2_local(:);
            
            % Remove NaN values from the local neighborhood
            valid_idx = ~isnan(U1_local) & ~isnan(U2_local);
            X_valid = X(valid_idx);
            Y_valid = Y(valid_idx);
            U1_valid = U1_local(valid_idx);
            U2_valid = U2_local(valid_idx);
            
            % Ensure there are enough valid points for affine fitting
            if numel(U1_valid) < 3 % At least 3 points 
                continue; % Skip smoothing for this point
            end
            
            % Construct design matrix for affine fit
            A = [X_valid, Y_valid, ones(length(X_valid), 1)];  
            
            % Solve least-squares for best-fit affine transformation
            beta_x = A \ U1_valid;  % Solve for x-displacement coefficients
            beta_y = A \ U2_valid;  % Solve for y-displacement coefficients
            
            % Compute smoothed displacement at (i, j)
            U1_smooth(i, j) = [j, i, 1] * beta_x;  % Smoothed U1
            U2_smooth(i, j) = [j, i, 1] * beta_y;  % Smoothed U2
        end
    end
    
    % Store results in output cell arrays
    u1_affine{c} = U1_smooth;
    u2_affine{c} = U2_smooth;
end

%% %%%%%%% CALCULATING DISPLACEMENT GRADIENTS %%%%%%%%%%% 

% Define grid spacing
dx = L / sz(2);
dy = H / sz(1);

% Initialize displacement gradient tensor components
H11 = cell(size(u1));
H12 = cell(size(u1));
H21 = cell(size(u1));
H22 = cell(size(u1));

for c = 1:num_cells

    % Empty node cells (bigger dimnesions)
    u1_node = NaN(N2+1, N1+1);
    u2_node = NaN(N2+1, N1+1);

    % Comptuing the nodes 
    for i = 2:N2
        for j = 2:N1
            if any(~isnan([u1_affine{c}(i,j), u1_affine{c}(i-1, j), u1_affine{c}(i, j-1), u1_affine{c}(i-1, j-1)]))
                u1_node(i, j) = mean([u1_affine{c}(i,j), u1_affine{c}(i-1, j), u1_affine{c}(i, j-1), u1_affine{c}(i-1, j-1)],'omitnan');
            end
            if any(~isnan([u2_affine{c}(i,j), u2_affine{c}(i-1, j), u2_affine{c}(i, j-1), u2_affine{c}(i-1, j-1)]))
                u2_node(i, j) = mean([u2_affine{c}(i,j), u2_affine{c}(i-1, j), u2_affine{c}(i, j-1), u2_affine{c}(i-1, j-1)],'omitnan');
            end
        end
    end

    % Numerical differentiation of the nodes
    H11{c}=NaN(N2,N1);
    for i = 1:N2
        for j = 1:N1
            if any(isnan([u1_node(i,j), u1_node(i+1, j), u1_node(i, j+1), u1_node(i+1, j+1)]))
                continue
            else
                H11{c}(i, j) = (0.5 * (u1_node(i, j+1) + u1_node(i+1, j+1))-0.5 * (u1_node(i,j) + u1_node(i+1, j))) / dx;
            end
        end    
    end
end


%Average H11 over the whole medium
H11_av=[];
for i=1:length(H11)
    H11_av(i)=mean(mean(H11{i},'omitnan'),'omitnan');
end

%Average only over x2
H11_x={};
for i=1:length(H11)
    H11_x{i}=mean(H11{i},'omitnan');
end

%% %%%%%%%%% CALCULATING THE EFFECTIVE STRESS %%%%%%%%%%%
points1=find(round(P_new.*10^-6,2)==0.15); %start point of PIV
p1=points1(1);
P_short=P_new(p1:end);

x1=dx/2:dx:L-dx/2; %Vector of x1

%Calculating the effective stress along x1
S={};
for i=1:length(P_short)
    S{i}=-(P_short(i)/L).*x1;
end

S_av=[];%Average effective stress over the whole meidum
for i=1:length(S)
    S_av(i)=mean(S{i});
end

%% %%%%%%%%% PLOTTING STRAiN-STRESS OVER THE WHOLE MEDIUM %%%%%%%%%%
point2=find(round(P_short.*10^-6,2)==0.02);
fig=figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.5, 0.6]); % Wider and shorter
hold on

place=find(round(S_av.*10^-6,6)==-0.212868); %Aproximated onset of localization (extracted manualy on the graph)
p=place(1);
classic_colors = colororder;
plot(S_av(1:point2(length(point2))).*10^-6,H11_av(1:point2(length(point2))).*10^2,'color',classic_colors(1,:),'Linewidth',2)
%5-sec window around the localization
x=[S_av(p+40).*10^-6 S_av(p-10).*10^-6 S_av(p-10).*10^-6 S_av(p+40).*10^-6];
y=[-1 -1 0 0 ];
patch(x,y,'k','FaceAlpha',0.1,'Linestyle','none','HandleVisibility','off')

ylabel('\epsilon [%]')
xlabel('\sigma` [MPa]')
ylim([-0.9 0])
xlim([-0.32 0])
set(gca,'FontSize',16)
xticks([-0.3 -0.15 0])
yticks([-0.9 -0.45 0])
set(gca, 'XDir', 'reverse')
grid on

%% %%%%%%%%%% PLOTTING STRIN AND STRESS ALONG THE MEDIUM %%%%%%%%%%
fig=figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.3, 0.7]); % Wider and shorter
t = tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');
p1=find(round(S_av.*10^-6,2)==-0.1);
p2=find(round(S_av.*10^-6,2)==-0.2);
p3=find(round(S_av.*10^-6,2)==-0.3);

places=[p1(1),p2(1),p3(1)];
classic_colors = colororder;
for i=places
    nexttile
    hold on
    plot(S{i}.*10^-6,H11_x{i}.*10^2,'color', classic_colors(1,:),'Linewidth',2);
    ylabel('\epsilon [%]')
    ylim([min(H11_x{i}.*10^2) 0])
    yticks([min(H11_x{i}.*10^2) 0])
    yticklabels({num2str(round(min(H11_x{i}).*10^2,1)), '0'})
    xlim([min(S{i}.*10^-6) 0])
    xticks([min(S{i}.*10^-6) mean(S{i}.*10^-6) 0])
    xticklabels({num2str(round(min(S{i}.*10^-6),1)), num2str(round(mean(S{i}.*10^-6),2)),'0'})
    grid on
    set(gca,'FontSize',16)
    set(gca, 'XDir', 'reverse')
end

nexttile(3)
xlabel('\sigma` [MPa]')

%% Saving the workspace
% close all
% save('E:\Compaction_localization_data\S1\Workspace_strain_stress_analyzed_S1.mat');