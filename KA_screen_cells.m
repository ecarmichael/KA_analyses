function KA_screen_cells(data_dir, out_dir)
%% KA_screen_cells: loop over the output figures from KA_screener in a dir
% and collect them in some tile plots for easier screening. 
%
%
%  Inputs:
%       data_dir: [string]  path to the intermediate files from
%       KA_screener

%% plot basics setup. 

c_ord = linspecer(5); 
ppf =4; 
R = 2;
C = 5; 

%% get the paths for all the rate plots. 

f_list = dir([data_dir filesep '**' filesep '*rate.fig']);

for iF = length(f_list) :-1:1
   if contains(f_list(iF).name, 'entry') 
       keep_idx(iF) = 0;
   else
       keep_idx(iF) = 1;
   end
end

f_list(~keep_idx) = []; 

for iF = 1:length(f_list)
    close all
    
    F1 = openfig([f_list(iF).folder filesep  f_list(iF).name]);
    subplot(1,3,3, 'parent', F1);
    H_r = gca;
    x_lims_r = xlim; 
    y_lims_r = ylim; 

    pause(1)
    
    TTID_idx = strfind(f_list(iF).name, 'rate');
    F2 = openfig([f_list(iF).folder filesep  f_list(iF).name(1:TTID_idx-1) 'all_zcore.fig']);
    subplot(2,1,1, 'parent',F2 );
    H_rast = gca; 
    y_lims = ylim; 
    
    subplot(2,1,2,'parent',F2 );

    yyaxis right
    H_gau_r = gca;
    velo = H_gau_r.Children;
    
    
    yyaxis left
    H_gau_l = gca;
            y_lims_gau = ylim; 


    f_main = figure(iF+100);
    
    S0 = subplot(3,2,1, 'parent', f_main); % thanks mvdm for this trick
    S1 = subplot(3,2,2, 'parent', f_main); % thanks mvdm for this trick
    S2 = subplot(3,2,3:4, 'parent', f_main); 
    S3 = subplot(3,2,5:6, 'parent', f_main); 

    fig1 = get(H_r,'children'); %get handle to all the children in the figure
    fig2 = get(H_rast,'children');
    fig4 = get(H_gau_l,'children');


    copyobj(fig1,S1); %copy children to new parent axes i.e. the subplot axes
    copyobj(fig2,S2);
    copyobj(fig4,S3);

    subplot(3,2,1,'parent', f_main)
            parts = strsplit(f_list(iF).folder, filesep);
    text(0, 0.6, {parts{end-1} ; f_list(iF).name(1:5)}, 'fontsize', 18);
    axis off
    
    
    subplot(3,2,2, 'parent', f_main); % thanks mvdm for this trick
        set(gca, 'YLim', y_lims_r, 'XLim',x_lims_r); 
axis off
%     xlim(H_r.XLim); ylim(H_r.YLim); axis off
    


   subplot(3,2,3:4, 'parent', f_main); 
    set(gca, 'YLim', y_lims); 
    
   subplot(3,2,5:6, 'parent', f_main); 
    set(gca, 'YLim', y_lims_gau); 
yyaxis right
plot(velo.XData, velo.YData)
        
        SetFigure([], gcf)
        
        set(gcf, 'position', [600    50   784   815]); 
        
        mkdir([out_dir filesep 'Screeners']);
        mkdir([out_dir filesep 'Screeners_fig']);
        saveas(gcf, [out_dir filesep 'Screeners' filesep parts{end-1} '_' f_list(iF).name(1:5) '.png'])
        saveas(gcf, [out_dir filesep 'Screeners_fig' filesep parts{end-1} '_' f_list(iF).name(1:5) '.png'])

end


