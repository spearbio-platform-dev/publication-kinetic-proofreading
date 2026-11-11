% Fig. S9. Requires the Symbolic Math Toolbox, fun_competitive_num.m and
% fun_PlotWithSimpleCI.m.
clear;

vfCA = logspace(-14,-6);
% vfCA = (0:10)*1e-9;
fSigma = 0.1;




%%
vfCB = [200e-12,   200e-12,   200e-12,  300e-12,   500e-12,    1e-9];
raCB = {'200 pM', '200 pM',   '200 pM', '300 pM',  '500 pM',  '1 nM'};
vfCE = [50e-12,    10e-12,     2e-12,   50e-12,    50e-12,    50e-12];
raCE = {'50 pM',  '10 pM',    '2 pM',   '50 pM',   '50 pM',   '50 pM'};
KD = 200e-12; sKD = '200 pM';
KDE = KD; %M

hF = figure('Units','points','Position',[10 10 1000 800]);
xiColors = get(gca,'colororder');
hT = tiledlayout(2,3);
vhA = gobjects(1,6);
hCounter = 0;

for ii = 1:6
    hCounter = hCounter + 1;
    vhA(hCounter) = nexttile;
    hA = gca;

    vfSig_True = fun_competitive_num(vfCA,[vfCB(ii), vfCE(ii), KD, KDE]);
    fSig_Bkgd = fun_competitive_num(0,[vfCB(ii), vfCE(ii), KD, KDE]);
    vfSig_Bkgd = vfSig_True - vfSig_True + fSig_Bkgd;
    iIdxGuess = find(vfSig_True*(1-fSigma) - fSig_Bkgd*(1+fSigma) > 0,1,'first');
    g = @(x)(fun_competitive_num(x,[vfCB(ii), vfCE(ii), KD, KDE])*(1-fSigma) - fSig_Bkgd*(1+fSigma));
    fLOD = fzero(g,vfCA(iIdxGuess));
    fun_PlotWithSimpleCI(vhA(hCounter),vfCA,vfSig_True,vfSig_True*(1+fSigma),vfSig_True*(1-fSigma),xiColors(2,:));
    fun_PlotWithSimpleCI(vhA(hCounter),vfCA,vfSig_Bkgd,vfSig_Bkgd*(1+fSigma),vfSig_Bkgd*(1-fSigma),xiColors(1,:));
    set(gca,'fontsize',12);    
    xlabel('$C_A$ (M)','Interpreter','latex');
    ylabel('$[E]_{eq}$ (M)','Interpreter','latex');
    title(...
        sprintf('$K_D$ = %s, $C_B$ = %s, $C_E$ = %s\n $\\Rightarrow$ LOD = %.1f pM',...
        sKD,raCB{ii},raCE{ii},fLOD*1e12),...
        'Interpreter','latex');
    set(gca,'Xlim',[min(vfCA) max(vfCA)]);
    set(gca,'Ylim',[1e-12 3e-10]);
    set(gca,'XScale','log','YScale','log');
    plot(vhA(hCounter),[fLOD, fLOD],vhA(hCounter).YLim,'--','color',[0.5 0.5 0.5]);
    legend({...
        '$S(C_A)$',...
        '$S(C_A)\cdot(1\pm\sigma)$',...
        '$S_{BG}$',...
        '$S_{BG}\cdot(1\pm\sigma)$',...
        'LOD'},...
        'Location','best',...
        'Interpreter','latex');
end
