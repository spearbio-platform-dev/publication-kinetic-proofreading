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

%%

vfKD = [200e-12, 1e-9 ];
raKD = {'200 pM','1 nM'};
vfCE = [...
    200e-12;
    50e-12;
    10e-12;
    2e-12];
raCB_FineGrain = cell(4,2);
raLOD = cell(4,2);
for i_KD = 1:length(vfKD)
    fKD = vfKD(i_KD);
    fKDE = fKD;
    for i_CE = 1:length(vfCE)
        fCE = vfCE(i_CE);
        vfCB_FineGrain = logspace(log10(fCE),log10(100*fCE));
        raCB_FineGrain{i_CE,i_KD} = vfCB_FineGrain;
        raLOD{i_CE,i_KD} = nan(1,50);
        for j = 1:length(vfCB_FineGrain)
            fprintf('i_CE = %u;j = %u\n',i_CE, j);
            fCB = vfCB_FineGrain(j);
            vfSig_True = fun_competitive_num(vfCA,[fCB, fCE, fKD, fKDE]);
            fSig_Bkgd = fun_competitive_num(0,[fCB, fCE, fKD, fKDE]);
            iIdxGuess = find(vfSig_True*(1-fSigma) - fSig_Bkgd*(1+fSigma) > 0,1,'first');
            if ~isempty(iIdxGuess)
                g = @(x)(fun_competitive_num(x,[fCB, fCE, fKD, fKDE])*(1-fSigma) - fSig_Bkgd*(1+fSigma));
                raLOD{i_CE,i_KD}(j) = fzero(g,vfCA(iIdxGuess));
            end
        end
    end
end
%%
for i_KD = 1:length(vfKD)    
    hCounter = hCounter + 1;
    vhA(hCounter) = nexttile;
    hold on;

    for i_CE = 1:length(vfCE)
        vfCB_FineGrain = raCB_FineGrain{i_CE,i_KD};
        vfLOD = raLOD{i_CE,i_KD};
        plot(vhA(hCounter),vfCB_FineGrain,vfLOD,'-o','LineWidth',2,'Color',xiColors(i_CE,:));
        plot([min(vfCB_FineGrain),max(vfCB_FineGrain)],[min(vfLOD),min(vfLOD)],'--','Color',xiColors(i_CE,:));
    end
    set(vhA(hCounter),'Fontsize',12);
    set(vhA(hCounter),'xlim',[min(cellfun(@(x)min(x),raCB_FineGrain)),max(cellfun(@(x)min(x),raCB_FineGrain))]);
    legend(vhA(hCounter),{'$C_E$ = 200 pM','','$C_E$ = 50 pM','','$C_E$ = 10 pM','','$C_E$ = 2 pM'},'Interpreter','latex');
    title(sprintf('$K_D$ = %s\n $\\Rightarrow$ min. LOD = %.1f pM',raKD{i_KD},min(min(xfLOD))*1e12),'Interpreter','latex');
    xlabel('$C_B$ (M)','Interpreter','latex');
    ylabel('LOD (M)','Interpreter','latex');
end



