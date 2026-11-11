% Figs. S8a (first figure) and S10 (second figure). The second figure also writes
% Allostery.mat (Fig. S10f data, last pass of the E_in/E_act loop).
% Requires the Symbolic Math Toolbox and fun_PlotWithSimpleCI.m.
clear;

%% Equilibrum Concentrations
syms A B Bs ABs KD K1 CA CB Sigma E0
eqns_Allostery = [...
    ABs/(A*Bs) == 1/KD,...
    A + ABs == CA,...
    Bs/B == K1,...
    B + Bs + ABs == CB];
SFn_Allostery = solve(eqns_Allostery,[A B Bs ABs]);


%% Compute LOD - B has some signal
eqns_LOD_Allo_whE0 = [eqns_Allostery,...
    (1/(1+K1)*CB*E0 +K1/(1+K1)*CB)*(1+Sigma) == (B*E0 + (Bs + ABs))*(1-Sigma)...
    ];
SFn_LOD_whE0 = solve(eqns_LOD_Allo_whE0,[A B Bs ABs CA]);
MFn_LOD_whE0 = matlabFunction(SFn_LOD_whE0.CA);

%% Additional parameters
fE0 = 0.01;
fCB = 0.2e-9;
fSigma = 0.1;
vfCA = logspace(-15,-6);
vfCA_Zero = vfCA - vfCA + 0;

MFn_B0  = matlabFunction(SFn_Allostery.B);
MFn_Bs  = matlabFunction(SFn_Allostery.Bs);
MFn_ABs = matlabFunction(SFn_Allostery.ABs);

vfKD = [200e-12, 1e-9, 5e-9];
raKD = {'200 pM','1 nM','5 nM'};
vfK1 = [1, 0.1, 0.01];

%% Get colors.

figure;
xiColors = get(gca,'colororder');
hA = gca;

fK1 = vfK1(1);
fKD = vfKD(1);

% Compute background
vfBs__Bkgd = MFn_Bs(vfCA_Zero,fCB,fK1,fKD);
vfBs__Bkgd = vfBs__Bkgd(vfBs__Bkgd >= 0)';
vfB0__Bkgd = MFn_B0(vfCA_Zero,fCB,fK1,fKD);
vfB0__Bkgd = vfB0__Bkgd(vfB0__Bkgd >= 0)';
vfSig_Bkgd = vfB0__Bkgd * fE0 + vfBs__Bkgd;


% Compute signal
vfB0__True = MFn_B0(vfCA,fCB,fK1,fKD);
vfB0__True = vfB0__True(vfB0__True >= 0)';
vfBs__True = MFn_Bs(vfCA,fCB,fK1,fKD);
vfBs__True = vfBs__True(vfBs__True >= 0)';
vfABs_True = MFn_ABs(vfCA,fCB,fK1,fKD);
vfABs_True = vfABs_True(vfABs_True >= 0 & vfABs_True <= vfCA & vfABs_True <= fCB)';
vfSig_True = vfB0__True * fE0 + vfBs__True + vfABs_True;

% Compute LOD
fLOD = MFn_LOD_whE0(fCB,fE0,fK1,fKD,fSigma);


fun_PlotWithSimpleCI(hA,vfCA,vfSig_True,vfSig_True*(1+fSigma),vfSig_True*(1-fSigma),xiColors(2,:));
fun_PlotWithSimpleCI(hA,vfCA,vfSig_Bkgd,vfSig_Bkgd*(1+fSigma),vfSig_Bkgd*(1-fSigma),xiColors(1,:));
xlabel('$C_A$ (M)','Interpreter','latex');
ylabel('Signal (eq. M)','Interpreter','latex');
set(gca,'Xlim',[1e-12 1e-8]);
set(gca,'Ylim',[6e-11 3e-10]);
set(gca,'XScale','log','YScale','log','fontsize',16);
plot(hA,[fLOD, fLOD],hA.YLim,'--','color',[0.5 0.5 0.5],'linewidth',2);
legend({...
    '$S(C_A)$',...
    '$S(C_A)\cdot(1\pm\sigma)$',...
    '$S_{BG}$',...
    '$S_{BG}\cdot(1\pm\sigma)$',...
    'LOD'},...
    'Location','best',...
    'Interpreter','latex');




%%
figure('Units','points','Position',[10 10 1000 800]);
tiledlayout(2,3);
fKD = vfKD(1);
sKD = raKD{1};

for i_K1 = 1:length(vfK1)
    nexttile;
    hA = gca;
    set(gca,'fontsize',12);
    fK1 = vfK1(i_K1);
    % Compute background
    vfBs__Bkgd = MFn_Bs(vfCA_Zero,fCB,fK1,fKD);
    vfBs__Bkgd = vfBs__Bkgd(vfBs__Bkgd >= 0)';
    vfB0__Bkgd = MFn_B0(vfCA_Zero,fCB,fK1,fKD);
    vfB0__Bkgd = vfB0__Bkgd(vfB0__Bkgd >= 0)';
    vfSig_Bkgd = vfB0__Bkgd * fE0 + vfBs__Bkgd;

    % Compute signal
    vfB0__True = MFn_B0(vfCA,fCB,fK1,fKD);
    vfB0__True = vfB0__True(vfB0__True >= 0)';
    vfBs__True = MFn_Bs(vfCA,fCB,fK1,fKD);
    vfBs__True = vfBs__True(vfBs__True >= 0)';
    vfABs_True = MFn_ABs(vfCA,fCB,fK1,fKD);
    vfABs_True = vfABs_True(vfABs_True >= 0 & vfABs_True <= vfCA & vfABs_True <= fCB)';
    vfSig_True = vfB0__True * fE0 + vfBs__True + vfABs_True;

    % Compute LOD
    fLOD = MFn_LOD_whE0(fCB,fE0,fK1,fKD,fSigma);

    fun_PlotWithSimpleCI(hA,vfCA,vfSig_True,vfSig_True*(1+fSigma),vfSig_True*(1-fSigma),xiColors(2,:));
    fun_PlotWithSimpleCI(hA,vfCA,vfSig_Bkgd,vfSig_Bkgd*(1+fSigma),vfSig_Bkgd*(1-fSigma),xiColors(1,:));
    xlabel('$C_A$ (M)','Interpreter','latex');
    ylabel('Signal (eq. M)','Interpreter','latex');
    title(sprintf('$K_D$ = %s, $K_1$ = %.2f\n $\\Rightarrow$ LOD = %.1f pM',sKD,fK1,fLOD*1e12),'Interpreter','latex');
    set(gca,'Xlim',[min(vfCA) max(vfCA)]);
    set(gca,'Ylim',[1e-12 3e-10]);
    set(gca,'XScale','log','YScale','log');
    plot(hA,[fLOD, fLOD],hA.YLim,'--','color',[0.5 0.5 0.5]);
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
for i_KD = 1:2
    fKD = vfKD(i_KD);
    sKD = raKD{i_KD};
    % LOD as a function of K1
    vfK1_fine = logspace(-3,0);
    vfLOD_fine = MFn_LOD_whE0(fCB,fE0,vfK1_fine,fKD,fSigma);
    
    nexttile;
    hold on;
    set(gca,'Fontsize',12);
    plot(vfK1_fine,vfLOD_fine,'LineWidth',2);
    plot([min(vfK1_fine),max(vfK1_fine)],[min(vfLOD_fine),min(vfLOD_fine)],'--','Color',[0.5 0.5 0.5]);
    title(sprintf('$K_D$ = %s\n $\\Rightarrow$ Min. LOD: %.1f pM\n',sKD,min(vfLOD_fine)*1e12),'Interpreter','latex');
    xlabel('$K_1$','Interpreter','latex');
    ylabel('LOD (M)','Interpreter','latex');
end

vfE0 = [0.01 0.001 0.0001];
nexttile;
hold on;
set(gca,'Fontsize',12);

for iE0 = 1:length(vfE0)
    fE0_var = vfE0(iE0);
    vfKD_fine = logspace(-13,-7);
    vfLOD_AllosteryFine = nan(size(vfKD_fine));
    for i_KD_fine = 1:length(vfKD_fine)
        fKD = vfKD_fine(i_KD_fine);
        % LOD as a function of K1
        vfK1_fine2 = logspace(-3,0);
        vfLOD_fine2 = MFn_LOD_whE0(fCB,fE0_var,vfK1_fine2,fKD,fSigma);
        vfLOD_AllosteryFine(i_KD_fine) = min(vfLOD_fine2);
    end
    plot(vfKD_fine,vfLOD_AllosteryFine,'LineWidth',2);
    
end

plot(vfKD_fine,vfKD_fine,'--');
plot(vfKD_fine,vfKD_fine.*2*fSigma,'--');
fill([1e-10 1e-8 1e-8 1e-10],[1e-7 1e-7 1e-14 1e-14],[0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor','none');

set(gca,'XScale','log','YScale','log',...
    'Xlim',[min(vfKD_fine) max(vfKD_fine)],'YLim',[1e-14 1e-7]);
set(gca,'XTick',logspace(-13, -7, 7),'YTick',logspace(-14,-7,8));
xlabel('$K_D$ (M)','Interpreter','latex');
ylabel('Min. LOD at optimal $K_1$ (M)','Interpreter','latex');
legend({'$E_{in}/E_{act}$ = 0.01','$E_{in}/E_{act}$ = 0.001','$E_{in}/E_{act}$ = 0.0001',...
    '$y$ = $x$','$y$ = 2$\sigma\cdot x$',sprintf('Typical $K_D$ for antibodies \n(100 pM to 10 nM)')},...
    'Interpreter','latex','Location','eastoutside');


sFileName = 'Allostery_v2.png';
saveas(gcf,sFileName);

save('Allostery.mat','vfKD_fine','vfLOD_AllosteryFine')
