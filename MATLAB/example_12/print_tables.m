% ----------------------------------------------------------------------- %
%  MM   MM   OOOOO   DDDDDD   EEEEEEE  LL         AAA    VV   VV   GGGGG  %
%  MM   MM  OOOOOOO  DDDDDDD  EEEEEEE  LL        AA AA   VV   VV  GG   GG %
%  MMM MMM  OO   OO  DD   DD  EE       LL        AA AA   VV   VV  GG   GG %
%  MM M MM  OO   OO  DD   DD  EEEEE    LL       AA   AA  VV   VV  GGGGGGG %
%  MM   MM  OO   OO  DD   DD  EEEEE    LL       AAAAAAA  VV   VV   GGGGGG %
%  MM   MM  OO   OO  DD   DD  EE       LL       AA   AA   VV VV        GG %
%  MM   MM  OOOOOOO  DDDDDDD  EEEEEEE  LLLLLLL  AA   AA    VVV     GGGGGG %
%  MM   MM   OOOOO   DDDDDD   EEEEEEE  LLLLLLL  AA   AA     V     GGGGGGG %
% ----------------------------------------------------------------------- %

% Look at results of evaluation period
table_res = table_res_val;

% Now print Table 8 using Latex language
table_content = 'heteros'; % 'homos_only'
row_table = {'$\text{QS}$','$\text{LS}$','$\text{SS}$','$\text{CRPS}$',...
    '$\text{ES}$','$R_\text{l}$','$C_\text{v}$','$C$','$W$',...
    '$\ell(\widehat{\bm{\Phi}}\vert \bm{\upomega})$','$\text{RMSE}$',...
    '$\text{NSE}$','\text{KGE}'};
% Which results [constant/nonconstant var + single/group treatment?]
switch table_content
    case 'homos'        % homoscedastic variance [group and single]
        % order     normal  logn    gnorm   tnorm  gamma   GEV
        idx_col = [  1 8  , 2 9   , 3 10  , 4 11  , 5 12  , 6 13 ];
    case 'heteros'      % heteroscedastic variance [group and single]
        % order     normal  logn    gnorm   tnorm   gamma    GEV
        idx_col = [ 15 22 , 16 23 , 17 24 , 18 25 , 19 26 , 20 27 ];
    case 'both_group'   % homos/heteros - with group variance treatment
        % order     normal  logn   gnorm    tnorm  gamma   GEV
        idx_col = [ 1 15 , 2  16 ,  3  17 , 4  18 , 5  19 , 6  20 ];
    case 'both_ind'     % homos/heteros - with group variance treatment
        idx_col = [ 8 22 , 9  23 ,  10 24 , 11 25 , 12 26 , 13 27 ];
end
nd = numel(idx_col);
% Now extract columns for Table 8 in paper
Table_8 = table_res(:,idx_col);
% Write in latex format
for i = 1:numel(row_table)
    str = cell(1,nd);
    for z = 1:nd
        if Table_8(i,z) > 0
            str{z} = strcat('\multicolumn{1}{r}{$',num2str(Table_8(i,z),'%5.3f'),'$}');
        else
            str{z} = strcat('\multicolumn{1}{r}{-$',num2str(abs(Table_8(i,z)),'%5.3f'),'$}');
        end
    end
    switch table_content
        case {'homos','heteros'}
            % order     normal  logn    gnorm   tnorm   gamma    GEV
            fprintf(['%s & %s & %s & %s & %s && %s & %s & %s & %s && %s' ...
                '& %s & %s & %s && %s & %s & %s & %s && %s & %s & %s &'...
                '%s && %s & %s & %s & %s \\\\'],...
                strcat(char(row_table(i))),str{1:end});
        otherwise
            % new order: lognormal, weibull (homos) normal, gen_normal, gamma (heteros)
            str_tab = strcat('%s');
            for ii = 1:nd
                str_tab = strcat(str_tab,' & %s');
            end
            str_tab = strcat(str_tab,'\\\\');
            fprintf(str_tab,strcat(char(row_table(i))),str{1:end});
    end
    fprintf('\n');
end

% Rank the distributions for each storing rule
format long
Rk = cell(1,5); ncol = numel(idx_col);
Rk_num = nan(ncol,5);
ii_col = [1 1 , 2 2 , 3 3 , 4 4 , 5 5 , 6 6];
for zz = 1:5
    [~,ii] = sort(-table_res(zz,idx_col));
    for u = 1:numel(ii)
        Rk{zz} = names_pdf(idx_col(ii));
    end
    % Now also numeric values
    Rk_num(1:ncol,zz) = idx_col(ii);
    Rk_num2(1:ncol,zz) = ii_col(ii);
end
% Sort based on likelihood
L_rank = cell(1); [~,ii] = sort(-table_res(10,:));
for u = 1:numel(ii)
    L_rank{1} = names_pdf(ii);
end

% Table 9 of paper
idx = (table_res(10,:)~=-inf);
R = corr(table_res(:,idx)'); R = [R(10,1:9) R(10,11:13)];
% Print Table 9 to screen
str = cell(1,11);
for z = 1:11
    if R(1,z) > 0
        str{z} = strcat('\multicolumn{1}{c}{$',num2str(R(1,z),'%5.3f'),'$}');
    else
        str{z} = strcat('\multicolumn{1}{c}{-$',num2str(abs(R(1,z)),'%5.3f'),'$}');
    end
end
str_tab = strcat('%s');
for ii = 1:11
    str_tab = strcat(str_tab,' & %s');
end
str_tab = strcat(str_tab,'\\\\');
fprintf(str_tab,str{1:end});
fprintf('\n');