%% %%%%%%%%% PREPARING THE DATA %%%%%%%%%%%
clc
clear
close all
datapath='E:\Compaction_localization_data\S3';
load(fullfile(datapath, 'Workspace_pre_localization_ref_S3_updated.mat'));
%load(fullfile(datapath, 'Workspace_pre_localization_ref_S3.mat'),'u_original','v_original');

% Orientation
u1={};
for i=1:length(u_original)
    u1{i}=flipud(fliplr(-u_original{i}));
end

u2={};
for i=1:length(v_original)
   u2{i}=flipud(fliplr(-v_original{i}));
end

%Dimensions
L=0.0286; %length in [m]
H=0.0088; %height in [m]
sz=size(u1{1});

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

%% Checking the smoothed displacement filed for t=0.01 [s]

figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.4]); 
% Create a tiled layout with no row gaps
tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact'); 

lim=[-3 3];
cm='jet';

nexttile
imagesc(u1_affine{1}.*10^6);
set(gca, 'YDir', 'normal')
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'FontSize', 16)
yticklabels({''})
xticklabels({''})
xticks([-1 sz(2)+1])
yticks([-1 sz(1)+1])
colormap(gca,cm)
caxis(lim) 

nexttile
imagesc(u1{1}.*10^6);
set(gca, 'YDir', 'normal')
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'FontSize', 16)
yticklabels({''})
xticklabels({''})
xticks([-1 sz(2)+1])
yticks([-1 sz(1)+1])
colormap(gca,cm)
caxis(lim) 



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

    H12{c}=NaN(N2,N1);
    for i = 1:N2
        for j = 1:N1
            if any(isnan([u1_node(i,j), u1_node(i+1, j), u1_node(i, j+1), u1_node(i+1, j+1)]))
                continue
            else
                H12{c}(i, j) = (0.5 * (u1_node(i+1, j) + u1_node(i+1, j+1))-0.5 * (u1_node(i,j) + u1_node(i, j+1))) / dy;
            end
        end    
    end

    H21{c}=NaN(N2,N1);
    for i = 1:N2
        for j = 1:N1
            if any(isnan([u2_node(i,j), u2_node(i+1, j), u2_node(i, j+1), u2_node(i+1, j+1)]))
                continue
            else
                H21{c}(i,j) = (0.5 * (u2_node(i, j+1) + u2_node(i+1, j+1))-0.5 * (u2_node(i,j) + u2_node(i+1, j))) / dx;
            end
        end    
    end

    H22{c}=NaN(N2,N1);
    for i = 1:N2
        for j = 1:N1
            if any(isnan([u2_node(i,j), u2_node(i+1, j), u2_node(i, j+1), u2_node(i+1, j+1)]))
                continue
            else
                H22{c}(i, j) = (0.5 * (u2_node(i+1, j) + u2_node(i+1, j+1))-0.5 * (u2_node(i,j) + u2_node(i, j+1))) / dy;
            end
        end    
    end
end
%% CALCULATING DISPLACEMENT GRADIENT TENSOR (H) ON BEST AFFINE FIT FIELDS

% Define grid spacing
dx = L / sz(2);
dy = H / sz(1);
Ny = sz(1);
Nx = sz(2);

% Initialize displacement gradient tensor components
H11 = {};
H12 = {};
H21 = {};
H22 = {};

for i = 1:length(u1_affine)

    % Transform into node deflections
    u_node = NaN(Ny+1, Nx+1);
    v_node = NaN(Ny+1, Nx+1);

    for ii = 2:Nx
        for jj = 2:Ny
            valid_vals = ~isnan([u1_affine{i}(jj, ii), u1_affine{i}(jj-1, ii), u1_affine{i}(jj, ii-1), u1_affine{i}(jj-1, ii-1)]);
            if sum(valid_vals) > 0
                u_node(jj, ii) = nanmean([u1_affine{i}(jj, ii), u1_affine{i}(jj-1, ii), u1_affine{i}(jj, ii-1), u1_affine{i}(jj-1, ii-1)]);
            end
            valid_vals = ~isnan([u2_affine{i}(jj, ii), u2_affine{i}(jj-1, ii), u2_affine{i}(jj, ii-1), u2_affine{i}(jj-1, ii-1)]);
            if sum(valid_vals) > 0
                v_node(jj, ii) = nanmean([u2_affine{i}(jj, ii), u2_affine{i}(jj-1, ii), u2_affine{i}(jj, ii-1), u2_affine{i}(jj-1, ii-1)]);
            end
        end
    end

    % Compute displacement gradient tensor components
    for j = 1:sz(2) % Columns
        for kk = 1:sz(1) % Rows
            if ~any(isnan([u_node(kk, j+1), u_node(kk+1, j+1), u_node(kk, j), u_node(kk+1, j)]))
                H11{i}(kk, j) = ( 0.5 * (u_node(kk, j+1) + u_node(kk+1, j+1)) - 0.5 * (u_node(kk, j) + u_node(kk+1, j)) ) / dx;
            else
                H11{i}(kk, j) = NaN;
            end
        end    
    end

    for j = 1:sz(2) % Columns
        for kk = 1:sz(1) % Rows
            if ~any(isnan([u_node(kk+1, j+1), u_node(kk+1, j), u_node(kk, j+1), u_node(kk, j)]))
                H12{i}(kk, j) = ( 0.5 * (u_node(kk+1, j+1) + u_node(kk+1, j)) - 0.5 * (u_node(kk, j+1) + u_node(kk, j)) ) / dy;
            else
                H12{i}(kk, j) = NaN;
            end
        end
    end

    for j = 1:sz(2) % Columns
        for kk = 1:sz(1) % Rows
            if ~any(isnan([v_node(kk, j+1), v_node(kk+1, j+1), v_node(kk, j), v_node(kk+1, j)]))
                H21{i}(kk, j) = ( 0.5 * (v_node(kk, j+1) + v_node(kk+1, j+1)) - 0.5 * (v_node(kk, j) + v_node(kk+1, j)) ) / dx;
            else
                H21{i}(kk, j) = NaN;
            end
        end    
    end

    for j = 1:sz(2) % Columns
        for kk = 1:sz(1) % Rows
            if ~any(isnan([v_node(kk+1, j+1), v_node(kk+1, j), v_node(kk, j+1), v_node(kk, j)]))
                H22{i}(kk, j) = ( 0.5 * (v_node(kk+1, j+1) + v_node(kk+1, j)) - 0.5 * (v_node(kk, j+1) + v_node(kk, j)) ) / dy;
            else
                H22{i}(kk, j) = NaN;
            end
        end
    end

end

%% %%%%%%% PLOTTING (NO H21) %%%%%

fig1=figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.7, 0.3]); 
% Create a tiled layout with no row gaps
t=tiledlayout(3,4, 'TileSpacing', 'compact', 'Padding', 'compact'); 

pos=[1 5 9 2 6 10 3 7 11 4 8 12];
n=1;
lim=[-0.15 0.15];
cm='jet';

for i=[1 2 10 200]
        nexttile(pos(n))
        imagesc(H11{i}.*10^2);
        set(gca, 'YDir', 'normal')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize', 16)
        yticklabels({''})
        xticklabels({''})
        colormap(gca,cm)
        caxis(lim) 


        nexttile(pos(n+1))
        imagesc(H12{i}.*10^2);
        set(gca, 'YDir', 'normal')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize', 16)
        yticklabels({''})
        xticklabels({''})
        colormap(gca,cm)
        caxis(lim) 

        nexttile(pos(n+2))
        imagesc(H22{i}.*10^2);
        set(gca, 'YDir', 'normal')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize', 16)
        yticklabels({''})
        xticklabels({''})
        colormap(gca,cm)
        caxis(lim) 

        n=n+3;

end

% Saving image
%exportgraphics(fig1, fullfile('E:\Compaction_localization_data\S3', 'Localization_S3.tiff'), 'Resolution', 600);

%% %%%%%%% PLOTTING (FULL) %%%%%

fig2=figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.4, 0.2]); 
% Create a tiled layout with no row gaps
t=tiledlayout(4,4, 'TileSpacing', 'compact', 'Padding', 'compact'); 

pos=[1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16];
n=1;
lim=[-0.15 0.15];
cm='jet';

for i=[1 2 10 200]

        nexttile(pos(n))
        imagesc(H11{i}.*10^2);
        set(gca, 'YDir', 'normal')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize', 16)
        yticklabels({''})
        xticklabels({''})
        colormap(gca,cm)
        caxis(lim) 


        nexttile(pos(n+1))
        imagesc(H12{i}.*10^2);
        set(gca, 'YDir', 'normal')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize', 16)
        yticklabels({''})
        xticklabels({''})
        colormap(gca,cm)
        caxis(lim) 

        nexttile(pos(n+2))
        imagesc(H21{i}.*10^2);
        set(gca, 'YDir', 'normal')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize', 16)
        yticklabels({''})
        xticklabels({''})
        colormap(gca,cm)
        caxis(lim) 

        nexttile(pos(n+3))
        imagesc(H22{i}.*10^2);
        set(gca, 'YDir', 'normal')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize', 16)
        yticklabels({''})
        xticklabels({''})
        colormap(gca,cm)
        caxis(lim) 

        n=n+4;
end

%Saving image
%exportgraphics(fig2, fullfile('E:\Compaction_localization_data\S3', 'Localization_S3_full.tiff'), 'Resolution', 600);
