# import numpy as np
# 
# from models.toolbox.dimensionless import Dimensionless
# from models.toolbox.geometry import Geometry
# from models.toolbox.hydraulics import Hydraulics
# from models.two_phase_flow.closures.stratified_closure import StratifiedClosure
# from models.two_phase_flow.closures.bubble_closure import BubbleClosure
# 
# 
# def sigmoidal(a,x1,x2):
#     return 1/(1+np.exp(-a*(x1-x2)))
# 
# class IntermittentClosure:
# 
#     @staticmethod
#     def model_closure(obj, t, qL, qG, alphaL, P):
# 
#         D = obj.parameters['D']()
#         g = obj.parameters['g']()
#         epw = obj.parameters['epw']()
#         rhoL = obj.parameters['rhoL']()
#         muL = obj.parameters['muL']()
#         rhoG = obj.parameters['rhoG']()
#         muG = obj.parameters['muG']()
#         tetha = obj.parameters['tetha']()
#         tension = obj.parameters['tension']()
# 
#         alphaG = 1-alphaL
# 
#         A = Geometry.area(D)
# 
#         vL = Hydraulics.velocity(qL, rhoL, alphaL*A)
#         vG = Hydraulics.velocity(qG, rhoG, alphaG*A)
# 
#         vLS = Hydraulics.velocity(qL, rhoL, A)
#         vGS = Hydraulics.velocity(qG, rhoG, A)
# 
#         vM = vGS + vLS
# 
#         ReL = Dimensionless.reynolds(D, vL, rhoL, muL)
#         ReLS = Dimensionless.reynolds(D, vLS, rhoL, muL)
#         ReM = Dimensionless.reynolds(D, vM, rhoL, muL)
#         
#         ## (0) Fator de transição: determina a curvatura da função que chaveia os
#         # regimes. Valores elevados promovem transição instantânea entre os regimes
#         # de escoamento.
#         fatorExp = 2e6;
#         
#         ## (1) Determinando parâmetros de entrada
#         tetha = parametros.tetha
#         D           = parametros.D;
#         epw        = parametros.epw;
#         hL          = parametros.h;
#         rhoL       = parametros.rhoL;
#         rhoG       = parametros.rhoG;
#         densL       = parametros.densL;
#         rhoG       = parametros.rhoG;
#         AL       = parametros.AL;
#         AG       = parametros.AG;
#         tension      = parametros.tension;
#         betah       = parametros.betah;
#         
#         # Calculando velocidades superficiais
#         vLS = gerav(qL,densL,AL+AG);
#         vGS = gerav(qG,rhoG,AL+AG);
#         
#         # Calculando velocidades reais
#         v_L = gerav(qL,densL,AL);
#         v_G = gerav(qG,rhoG,AG);
#         
#         # Calculando valores médios
#         rhoM = (1-alphaL).*densL+alphaL.*rhoG;
#         v_M = gerav(qL+qG,rhoM,AL+AG);
#         
#         # Calculando fração de área de líquido
#         alphaLL = 1-alphaL;
#         
#         # Gerando flag que identifica se fluidos estão ou não em sentido opostos
#         flag_SentidoOposto = (sign(v_L.*v_G)+1)/2;
#         
#         # Fatores de atrito de fanning segundo velocidades superficiais
#         [fS_L,ReS_L]=fator(vLS,densL,rhoL,D,epw);
#         [fS_G,ReS_G]=fator(vGS,rhoG,rhoG,D,epw);
#         
#         # Fatores de atrito de fanning segundo velocidades reais
#         [f_L,Re_L]=fator(v_L,densL,rhoL,D,epw);
#         [f_G,Re_G]=fator(v_G,rhoG,rhoG,D,epw);
#         
#         # Fator de atrito médio
#         f_M = (1-alphaL).*f_L+alphaL.*f_G;
#         
#         ## (2) Calculando parametros
#         
#         # Calculando grupos adimensionais (FG e FL)
#         [FL, FG] = geraFroude(densL, rhoG, D, vLS, vGS);
#         
#         # Calculando FG/(cos(teta))^0.5
#         FGast = FG./((cos(tetha)).^0.5);
#         
#         #    Calculando NL, K e T
#         NL = (fS_L.*FL.^2./cos(tetha)).^0.5;
#         K  = FG.*(ReS_L./cos(tetha)).^0.5;
#         T  = (2.*fS_L.*FL.^2./cos(tetha)).^0.5;
#         
#         EOD = D.*((densL-rhoG).*9.81./tension).^0.5;
#         
#         # Calculando fatores criticos
#         FG_c = (75./(EOD.^2.*((1+ (EOD./75).^2./f_L).^0.5-1))).^0.5;
#         FGast_c      = (((1-hL./D).^2).*(alphaL.^3).*(np.pi.*D./4)./(D.*sin(betah./2))).^0.5;
#         K_c          = 2.*alphaLL.^0.5.*alphaL./(0.01.^0.5);
#         FL_c         = 1.5.*(hL./D).^0.5.*alphaLL;
#         NL_c         = 1.0.*(1-hL./D).^0.5.*alphaLL./(f_L./fS_L).^0.5;
#         T_c          = (2.*np.pi.*(alphaLL.^2.*alphaL)./sin(betah./2)./(f_L./fS_L)).^0.5;
#         
#         
#         ## (3) Assumindo-se estratificado, determinando relação hL/D (VÁ PARA ITEM
#         
#         ## (4) Testando instabilidade Kelvin-Helmholtz
#         #    Para o caso de FG/(cos(teta))^0.5 maior que crítico -> (VÁ ITEM 8)
#         #    Para o caso de FG/(cos(teta))^0.5 menor ou igual ao crítico -> (VÁ PARA ITEM 5)
#         
#         [flag_naoest1] = sigmoidal(fatorExp,FGast,FGast_c);# 0 - se menor e 1 se maior
#         
#         ## (5) Verificando sentido
#         #    Fluido subindo (VÁ PARA ITEM 6)
#         #    Fluido descendo (VÁ PARA ITEM 5.1)
#         [flag_sobe]=(sign(cos(tetha).*(qL+qG))+1)./2; # 0 se desce e 1 - se sobe
#         flag_desce = 1-flag_sobe;
#         
#         ## (5.1) Testando a trajetória do filme líquido
#         #    Para o caso de NL maior que crítico -> bubble, intermitente, anular ou anular disperso (VÁ ITEM 8)
#         #    Para o caso de NL menor ou igual ao crítico -> estratificado smooth ou wavy (VÁ ITEM 7)
#         [flag_naoest2] = sigmoidal(fatorExp,NL,NL_c).*flag_desce;# 0 - se menor e 1 se maior
#         
#         ## (6) Testando iterações de Jeffreus (vento-onda) para escoamento
#         # estratificado comondas
#         #    Para o caso de K maior que crítico -> estratificado wavy (FIM)
#         #    Para o caso de K menor ou igual a crítico -> estratificado smooth (FIM)
#         [flag_wavy1] = sigmoidal(fatorExp,K,K_c);# 0 - se menor e 1 se maior
#         
#         ## (7) Testando ondas de gravidade
#         #    Para o caso de FL maior que crítico e fluido descendo -> estratificado wavy (FIM)
#         #    Para o caso contrário -> estratificao smooth (FIM)
#         [flag_wavy2] = sigmoidal(fatorExp,FL,FL_c);# 0 - se menor e 1 se maior
#         
#         ## (8) (Vem do 5.1) Assumindo-se anular#####################################
#         
#         ## (9) Testando bridging do líquido acima de abaixo de alphaLL = 0.35
#         
#         #    Se alphaLL>0.35 -> bubble ou intermitente (VÁ PARA ITEM 11)
#         flag_naoanular = sigmoidal(fatorExp,alphaLL,0.35);
#         
#         #    Se teta=90º -> bubble, intermitente, anular ou anular disperso (VÁ PARA ITEM 10)
#         flag_vertical = (abs(tetha)-90)<(10/180*np.pi); 
#         
#         #    Caso contrário -> anular ou anular disperso (FIM)
#         
#         ## (10) (Vem do 9) Testando estabilidade do filme para escoamento vertical
#         #    Para o caso de FG menor que crítico -> bubble ou intermitente (VÁ PARA ITEM 11)
#         #    Para o caso de FG maior ou igual a crítico -> anular ou anular disperso (FIM)
#         flag_naohorizontal = (abs(tetha))>(10/180*np.pi); 
#         flag_anular1 = sigmoidal(fatorExp,FG,FG_c)+flag_naohorizontal;# 0 - se menor e 1 se maior
#         flag_anular1(find(flag_anular1==2))=1; 
#         
#         ## (11) Testando bubble ou intermitente###################################
#         
#         #    Para o caso de |teta|<10º (VÁ PARA ITEM 12)
#         flag_horizontal = abs(tetha)<(10/180*np.pi);
#         #    Se alphaG>0.52 -> intermitente (FIM)
#         flag_intermitente1 = sigmoidal(fatorExp,alphaL,0.52);
#         #    Se alphaG<0.25 -> dispersed bubble (FIM)
#         flag_bubble1 = sigmoidal(fatorExp,0.25,alphaL);
#         #    Se 0.25<alphaG<0.52 -> determinar tamanhos das bolhas (ver equações)
#         dmax = (0.725+4.15.*(vGS./v_L).^0.5).*(tension./densL).^(3/5).*(2.*f_M.*v_M.^3./D).^(-2/5);
#         ddef = 2.*(0.4.*tension./(densL-rhoG)./9.81).^0.5;
#         dmigr = 3./8.*densL./(densL-rhoG).*f_M.*v_M.^2./9.81./cos(tetha);
#         flag_checkbolhas = (1-flag_bubble1).*(1-flag_intermitente1);
#         #    Se dmax<ddef e dmax< dmigr -> dispersed bubble (FIM)
#         flag_bubble2 = sigmoidal(fatorExp,ddef,dmax).*sigmoidal(fatorExp,dmigr,dmax).*flag_checkbolhas;
#         #    Se dmax>ddef ou dmax< dmigr -> intermitente (FIM)
#         
#         
#         ## (12) Testando empuxo vs. flutuações turbulentas
#         
#         #    Para o caso de T > crítico -> bubble (FIM)
#         #    Para o caso de T < crítico -> intermitente (FIM)
#         [flag_bubble3] = sigmoidal(fatorExp,T,T_c);# 0 - se menor e 1 se maior
#         
#         # Chaveadores principais
#         fl.Est = (1-flag_naoest1).*flag_sobe + (1-flag_naoest2).*(1-flag_sobe);
#         
#         fl.EsW = fl.Est.*(flag_wavy1.*flag_sobe+flag_wavy2.*(1-flag_sobe));
#         
#         fl.EsS = fl.Est.*(1-fl.EsW);
#         
#         fl.Anu = (1-fl.Est).*...
#                  (1-flag_naoanular).*...
#                  (flag_vertical.*flag_anular1);
#         
#         fl.Bub = (1-fl.Est).*...
#                  (1-fl.Anu).*...
#                  (flag_bubble3.*flag_horizontal+(flag_bubble1+flag_bubble2).*(1-flag_horizontal));
#         
#         fl.Int = (1-fl.Est).*...
#                  (1-fl.Anu).*...
#                  (1-fl.Bub);
# 
#  # Armazendo Fatores
# fatores.FL = FL;
# fatores.FG = FG;
# fatores.NL = NL;
# fatores.K  = K;
# fatores.T  = T;
# fatores.FGcosteta = FGast;
# 
#  # Armazendo Fatores Críticos
# fatoresC.FL = FL_c;
# fatoresC.NL = NL_c;
# fatoresC.K  = K_c;
# fatoresC.T  = T_c;
# 
# # Gerando marcadores de regime
# reg1 = fl.Anu.*1+fl.Bub.*2+fl.Int.*3+fl.Est.*4;
# reg2 = fl.Anu.*1+fl.Bub.*2+fl.Int.*3+fl.EsS.*4+fl.EsW.*5;
#         
#         
# 
#         return gammaL, gammaG, gammaI, dPL, dPG, dPI
