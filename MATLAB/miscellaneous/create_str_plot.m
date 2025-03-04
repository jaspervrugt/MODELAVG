function [str_table,str_plot] = create_str_plot(method,options,K)
% Generates string for output Tables and postprocessor figures

% First print only the beta values
str = cell(1,K);
for i = 1:K
    str{i} = ['beta_{',num2str(i),'}'];
end

% Now for BMA method expand the vector
if strcmp(method,'bma')
    switch options.VAR
        case '1'
            str{K+1} = 'sigma';
        case '2'
            for i = 1:K
                str{K+i} = ['sigma_{',num2str(i),'}'];
            end
        case '3'
            str{K+1} = 'c';
        case '4'
            for i = 1:K
                str{K+i} = ['c_{',num2str(i),'}'];
            end
    end
    K2 = numel(str);
    if sum(strcmp(options.PDF,{'gen_normal','gev','gpareto'}))
        if strcmp(options.PDF,'gen_normal')
            switch options.TAU
                case '1'
                    str{K2+1} = 'tau';
                case '2'
                    for i = 1:K
                        str{K2+i} = ['tau_{',num2str(i),'}'];
                    end
            end
        end
        if strcmp(options.PDF,'gev')
            switch options.TAU
                case '1'
                    str{K2+1} = 'xi';
                case '2'
                    for i = 1:K
                        str{K2+i} = ['xi_{',num2str(i),'}'];
                    end
            end
        end
        if strcmp(options.PDF,'gpareto')
            switch options.TAU
                case '1'
                    str{K2+1} = 'k';
                case '2'
                    for i = 1:K
                        str{K2+i} = ['k_{',num2str(i),'}'];
                    end
            end
        end
    end
end

% Table does not need latex print
str_table = str;
% Now make up each cell for latex printing
% beta --> $\beta$, c --> $c$, xi --> $xi$, k --> $k$, etc. 
str_plot = strrep(str,'beta','\beta');
str_plot = strrep(str_plot,'sigma','\sigma');
str_plot = strrep(str_plot,'tau','\tau');
str_plot = strrep(str_plot,'xi','\xi');
for i = 1:numel(str_plot)
    str_plot{i} = ['$' str_plot{i} '$'];
end

% % str = strcat('{','''beta_{1}'',');
% % for i = 2 : K-1
% %     str = strcat(str,'''beta_{',num2str(i),'}'',');
% % end
% % 
% % % Check methods
% % if any(strcmp(method,{'ewa','bga','aica','bica','gra','mma','mma-s'})) % sum(strncmp(method,{'ewa','bga','aica','bica','gra','mma','mma-s'},inf))
% %     % Create string for printing of weights to screen
% %     str = strcat(str,'''beta_{',num2str(K),'}''}');
% % end
% % 
% % % Now check out BMA
% % if strcmp(method,'bma')
% %     str = strcat(str,'''beta_{',num2str(K),'}'',');
% %     switch options.VAR
% %         case '1'
% %             str = strcat(str,'''sigma''');
% %         case '2'
% %             for i = 1:K-1
% %                 str = strcat(str,'''sigma_{',num2str(i),'}'',');
% %             end
% %             str = strcat(str,'''sigma_{',num2str(K),'}''');
% %         case '3'
% %             str = strcat(str,'''c''');
% %         case '4'
% %             for i = 1:K-1
% %                 str = strcat(str,'''c_{',num2str(i),'}'',');
% %             end
% %             str = strcat(str,'''c_{',num2str(K),'}''');
% %     end
% %     if sum(strcmp(options.PDF,{'gen_normal','gev','gpareto'}))
% %         if strcmp(options.PDF,'gen_normal')
% %             switch options.TAU
% %                 case '1'
% %                     str = strcat(str,',''tau','''}');
% %                 case '2'
% %                     for i = 1:K-1
% %                         str = strcat(str,',''tau_{',num2str(i),'}''');
% %                     end
% %                     str = strcat(str,',''tau_{',num2str(K),'}''}');
% %             end
% %         end
% %         if strcmp(options.PDF,'gev')
% %             switch options.TAU
% %                 case '1'
% %                     str = strcat(str,',''xi','''}');
% %                 case '2'
% %                     for i = 1:K-1
% %                         str = strcat(str,',''xi_{',num2str(i),'}''');
% %                     end
% %                     str = strcat(str,',''xi_{',num2str(K),'}''}');
% %             end
% %         end
% %         if strcmp(options.PDF,'gpareto')
% %             switch options.TAU
% %                 case '1'
% %                     str = strcat(str,',''k','''}');
% %                 case '2'
% %                     for i = 1:K-1
% %                         str = strcat(str,',''k_{',num2str(i),'}''');
% %                     end
% %                     str = strcat(str,',''k_{',num2str(K),'}''}');
% %             end
% %         end
% %     else
% %         str = strcat(str,'}');
% %     end
% % end
% % 
% % % Now we store this string for print to screen later
% % str_plot = str; str_plot = strrep(str_plot,'beta','\beta');
% % str_plot = strrep(str_plot,'sigma','\sigma');
% % str_plot = strrep(str_plot,'tau','\tau');
% % str_plot = strrep(str_plot,'xi','\xi');
% % str_plot = strrep(str_plot,'k','$k');
% % str_plot = strrep(str_plot,'\be','$\be');
% % str_plot = strrep(str_plot,'\sig','$\sig');
% % str_plot = strrep(str_plot,'\ta','$\ta');
% % str_plot = strrep(str_plot,'\xi','$\xi');
% % for i = 1 : K
% %     str_plot = strrep(str_plot,strcat(num2str(i),'}'),strcat(num2str(i),'}$'));
% % end
% % 
% % if isfield(options,'VAR')
% %     switch options.VAR
% %         case '1'
% %             str_plot = strrep(str_plot,'ma','ma$');
% %         case '2'
% %             % nothing needed
% %         case '3'
% %             str_plot = strrep(str_plot,'c','$c$');
% %         case '4'
% %             str_plot = strrep(str_plot,'c','$c');
% %     end
% % end
% % if isfield(options,'TAU')
% %     if strcmp(options.PDF,'gen_normal')
% %         switch options.TAU
% %             case '1'
% %                 str_plot = strrep(str_plot,'au','au$');
% %             case '2'
% %                 % nothing needed
% %         end
% %     end
% %     if strcmp(options.PDF,'gev')
% %         switch options.TAU
% %             case '1'
% %                 str_plot = strrep(str_plot,'xi','xi$');
% %             case '2'
% %                 % nothing needed
% %         end
% %     end
% %     if strcmp(options.PDF,'gpareto')
% %         switch options.TAU
% %             case '1'
% %                 str_plot = strrep(str_plot,'k','k$');
% %             case '2'
% %                 % nothing needed
% %         end
% %     end
% % end
% % str_plot = eval(str_plot);
% % % Now evaluate regular print string
% % str_table = eval(str);