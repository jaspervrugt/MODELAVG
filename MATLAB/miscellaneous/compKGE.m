function KGE = compKGE(Yobs,Ysim)
% This function computes the Kling-Gupta efficiency

method = 1;
switch method 
    case 1 % Alternative form II
        m_X = mean(Ysim); s_X = std(Ysim);  
        m_Y = mean(Yobs); s_Y = std(Yobs); 
        alfa = m_X/m_Y; beta = s_X/s_Y; 
        m_XY = mean(Ysim.*Yobs);
        r = (m_XY - m_X*m_Y)/(s_X*s_Y);
        term = 3 + alfa^2 + beta^2 + r^2 - 2 * (alfa + beta + r);
        KGE = 1 - sqrt(term);   
    case 2 % Alternative form I
        m_X = mean(Ysim); m2_X = m_X^2; s_X = std(Ysim); s2_X = var(Ysim);
        m_Y = mean(Yobs); m2_Y = m_Y^2; s_Y = std(Yobs); s2_Y = var(Yobs);
        m_XY = mean(Ysim.*Yobs); m2_XY = m_XY^2;
        term = 3 + m2_X/m2_Y - 2*m_X/m_Y + s2_X/s2_Y - 2*s_X/s_Y ...
            + (m2_XY + m2_X*m2_Y - 2*m_X*m_Y*m_XY) / (s2_X*s2_Y) ...
            - 2*(m_XY - m_X*m_Y)/(s_X*s_Y);
        % Compute KG efficiency
        KGE = 1 - sqrt(term);
    case 0 % Old form
        a = Yobs; b = Ysim;
        s_meas = std(a); m_meas = mean(a); s_sim = std(b);
        m_sim = mean(b); alpha = s_sim/s_meas; gamma = m_sim/m_meas;
        %r = corrcoef(b,a); r = r(1,2);
        r = ( mean(Ysim.*Yobs) - mean(Ysim)*mean(Yobs) ) / (std(Ysim)*std(Yobs));
        % Compute KG efficiency
        KGE = 1 - sqrt((r-1)^2 + (alpha-1)^2 + (gamma-1)^2);
end

% % a = Yobs; b = Ysim;
% % % compute components
% % Ea = mean(a); E2a = Ea^2; stda = std(a); Vara = var(a);
% % Eb = mean(b); E2b = Eb^2; stdb = std(b); Varb = var(b);
% % Eab = mean(a.*b); E2ab = Eab^2;
% second term
% term2a = (E2a + E2b - 2*Ea*Eb)/E2a;
% term2b = (Vara + Varb - 2*stda*stdb)/Vara;
% term2c = (Eab - Ea*Eb - stda*stdb)^2/(Vara*Varb);
% New try
% % term2a = (E2a + E2b - 2*Ea*Eb)/E2a;
% % term2b = (Vara + Varb - 2*stda*stdb)/Vara;
% % term2c = (E2ab + E2a*E2b + Vara*Varb -2*Eab*Ea*Eb - 2*Eab*stda*stdb + 2*Ea*Eb*stda*stdb)/(Vara*Varb);
% % term2bc = ( (Vara*Varb + Varb*Varb - 2*stda*stdb*Varb) + (E2ab + E2a*E2b + Vara*Varb -2*Eab*Ea*Eb - 2*Eab*stda*stdb + 2*Ea*Eb*stda*stdb) ) / (Vara*Varb)
% % term2abc = ( (E2a*Vara*Varb + E2b*Vara*Varb - 2*Ea*Eb*Vara*Varb) + ...
% %     (E2a*Vara*Varb + E2a*Varb*Varb - 2*E2a*stda*stdb*Varb) + ...
% %     (E2a*E2ab + E2a*E2a*E2b + E2a*Vara*Varb - 2*E2a*Eab*Ea*Eb - 2*E2a*Eab*stda*stdb + 2*E2a*Ea*Eb*stda*stdb) ) / (E2a*Vara*Varb);
% % term2abc = ( ( 3*E2a*Vara*Varb + E2b*Vara*Varb - 2*Ea*Eb*Vara*Varb) + ...
% %     (E2a*Varb*Varb - 2*E2a*stda*stdb*Varb) + ...
% %     (E2a*E2ab + E2a*E2a*E2b - 2*E2a*Eab*Ea*Eb - 2*E2a*Eab*stda*stdb + ...
% %     2*E2a*Ea*Eb*stda*stdb) ) / (E2a*Vara*Varb);
% % KGE = 1 - sqrt(term2abc);
% compute components
% % m_X = mean(Ysim); m2_X = m_X^2; s_X = std(Ysim); s2_X = var(Ysim);
% % m_Y = mean(Yobs); m2_Y = m_Y^2; s_Y = std(Yobs); s2_Y = var(Yobs);
% % m_XY = mean(a.*b); m2_XY = m_XY^2;

% % term2 = (m2_X*s2_X*s2_Y + 3*m2_Y*s2_X*s2_Y - 2*m_X*m_Y*s2_X*s2_Y ...
% %     + m2_Y*s2_X^2 - 2*s_X*s_Y*m2_Y*s2_X + m2_Y*m2_XY + m2_X*m2_Y^2 ...
% %     - 2*m_XY*m_X*m_Y^3 - 2*m_XY*m2_Y*s_X*s_Y + 2*m_Y^3*m_X*s_X*s_Y) ...
% %     / (m2_Y*s2_X*s2_Y); 
% % KGE = 1 - sqrt(term2);

% Alternative
% % term3 = m2_X/m2_Y + 3 - 2*m_X/m_Y + s2_X/s2_Y - 2*s_X/s_Y ...
% %     + m2_XY/(s2_X*s2_Y) + m2_X*m2_Y/(s2_X*s2_Y) ...
% %     - 2*m_X*m_Y*m_XY/(s2_X*s2_Y) - 2*m_XY/(s_X*s_Y) + 2*m_X*m_Y/(s_X*s_Y);
% % KGE = 1 - sqrt(term3);


