import numpy as np
import math, sys, os, time
count1 = time.perf_counter()
D = 20/1000; Dcu = 10/1000; Dlng = 10/1000
L = 1000/1000; Ma = 0.01220158; Vlng = 0.0
Tinlet = 80; Pinlet = 2    #atm
Iop = 1810; m = 50
n = 1000 
dx = L / n
totalT = int(3000)
Qp = np.zeros(n+1)

Pholng = 486.55; Dvislng = 0.000333; Lbdlng = 0.234; Cplng = 3077
Relng = Vlng * Dlng * Pholng / Dvislng; Prlng = Dvislng * Cplng / Lbdlng
hlng = 0
hLNG = 0.023 * Relng ** 0.8 * Prlng ** 0.4 * (Lbdlng / Dlng)
print(hlng,hLNG)

def cal_N2(P):
    global Ts, HL, latent, Srften, Dvisg, Dvisl, Phog, Phol, Lbdg, Lbdl, Cpg, Cpl
    if P <= 5.5 :   

        ap = 67.3176927096
        bp = 12.546000654
        cp = -2.9470223048
        dp = 0.4260771126
        ep = -0.02499361
        Ts = ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4

        ap = 1.9994652761
        bp = 0.0438484777
        cp = -0.0027600585
        dp = 0.0005022281
        ep = -0.0000264141
        Cpl = (ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4) * 1000

        ap = 1.0586029782
        bp = 0.0703775442
        cp = -0.006729199
        dp = 0.001077922
        ep = -0.0000584154
        Cpg = (ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4) * 1000

        ap = 164.6068761733
        bp = -24.8371177839
        cp = 5.8552706446
        dp = -0.848648198
        ep = 0.0498870955
        Lbdl = (ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4) * 0.001

        ap = 6.0345504887
        bp = 1.4013177848
        cp = -0.2949515939
        dp = 0.0436671607
        ep = -0.0025578509
        Lbdg = (ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4) * 0.001

        ap = 851.0809084429
        bp = -55.3392384078
        cp = 12.2036189672
        dp = -1.7610591039
        ep = 0.1030202265
        Phol = ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4

        ap = 0.2655932558
        bp = 4.429713006
        cp = -0.1586212831
        dp = 0.0239289639
        ep = -0.0012683308
        Phog = ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4

        ap = 225.8930317449
        bp = -88.3807648444
        cp = 26.9063106471
        dp = -4.2114538926
        ep = 0.256562992
        Dvisl = (ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4) * 0.000001

        ap = 4.6720598935
        bp = 0.9523090776
        cp = -0.2127034238
        dp = 0.0308504303
        ep = -0.0018075802
        Dvisg = (ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4) * 0.000001

        ap = 11.1868409505
        bp = -2.9174328461
        cp = 0.709817647
        dp = -0.1030590387
        ep = 0.0060574631
        Srften = ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4

        ap = 211.7309943639
        bp = -15.1021768565
        cp = 3.026715853
        dp = -0.4359143835
        ep = 0.0254535023
        latent = (ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4) * 1000

        ap = -152.6019039
        bp = 23.98134931
        cp = -5.581298165
        dp = 0.809115066
        ep = -0.047502666
        HL = (ap + bp * P + cp * P ** 2 + dp * P ** 3 + ep * P ** 4) * 1000 

def cal_Ps(T):
    if T <= 100:
    #78-108
        a = -16.672172091903
        b = 0.742515139959629
        c = -1.03339686599366E-02
        d = 2.7791531470053E-05
        e = 2.57515024734318E-07
    else:
    #100-cri
        a = 465.887016661309
        b = -17.0456484631514
        c = 0.235441831473843
        d = -1.48084879916672E-03
        e = 3.7289710165672E-06
    Ps = a + b * T + c * T ** 2 + d * T ** 3 + e * T ** 4
    return Ps

def cal_Cu(T):
    global Lbdcu, Cpcu, Phocu, Rcu
    Lbdcu = 400
    Cpcu = 380    
    Phocu = 9000
    Rcu = 1*10**(-8)

def cal_TF():
    global TONB, TFL, qCHF
    P0 = Pinlet
    Pho1Cr0 = 8462.600562; Lbd1Cr0 = 11.61052231; Cp1Cr0 = 371.8244787
    Ts0 = 67.3176927096 + 12.546000654 * P0 - 2.9470223048 * P0 ** 2 + 0.4260771126 * P0 ** 3 - 0.02499361 * P0 ** 4
    Cpl0 = (1.9994652761 + 0.0438484777 * P0 - 0.0027600585 * P0 ** 2 + 0.0005022281 * P0 ** 3 - 0.0000264141 * P0 ** 4) * 1000
    Cpg0 = (1.0586029782 + 0.0703775442 * P0 - 0.006729199 * P0 ** 2 + 0.001077922 * P0 ** 3 - 0.0000584154 * P0 ** 4) * 1000
    Lbdl0 = (164.6068761733 - 24.8371177839 * P0 + 5.8552706446 * P0 ** 2 - 0.848648198 * P0 ** 3 + 0.0498870955 * P0 ** 4) * 0.001
    Lbdg0 = (6.0345504887 + 1.4013177848 * P0 - 0.2949515939 * P0 ** 2 + 0.0436671607 * P0 ** 3 - 0.0025578509 * P0 ** 4) * 0.001
    Phol0 = 851.0809084429 - 55.3392384078 * P0 + 12.2036189672 * P0 ** 2 - 1.7610591039 * P0 ** 3 + 0.1030202265 * P0 ** 4
    Phog0 = 0.2655932558 + 4.429713006 * P0 - 0.1586212831 * P0 ** 2 + 0.0239289639 * P0 ** 3 - 0.0012683308 * P0 ** 4
    Dvisl0 = (225.8930317449 - 88.3807648444 * P0 + 26.9063106471 * P0 ** 2 - 4.2114538926 * P0 ** 3 + 0.256562992 * P0 ** 4) * 0.000001
    Dvisg0 = (4.6720598935 + 0.9523090776 * P0 - 0.2127034238 * P0 ** 2 + 0.0308504303 * P0 ** 3 - 0.0018075802 * P0 ** 4) * 0.000001
    Srften0 = (11.1868409505 - 2.9174328461 * P0 + 0.709817647 * P0 ** 2 - 0.1030590387 * P0 ** 3 + 0.0060574631 * P0 ** 4) / 1000
    latent = (211.7309943639 - 15.1021768565 * P0 + 3.026715853 * P0 ** 2 - 0.4359143835 * P0 ** 3 + 0.0254535023 * P0 ** 4) * 1000
    #TONB = Ts0 + 2 + 0.0071 * P0
    TFL = Ts0 + (126.1 - Ts0) / (88 * Dvisl0 ** 0.5) * (1 + ((Lbdl0 * Phol0 * Cpl0) / (Lbd1Cr0 * Pho1Cr0 * Cp1Cr0)) ** 0.5)
    Wel = ((Ma / (3.14 * D * D / 4)) ** 2 * L) / (Phol0 * Srften0)
    #qCHF = 0.486071 * (Ma / (3.14 * D * D / 4)) * latent * (Phog0 / Phol0) ** 0.6353 * Wel ** (-0.4987) * (L / D) ** 0.3112
    qCHF = 100000

def cal_l(ii):
    global Ma, Ml, Mg, TONB, TFL, qCHF
    V_temp = Ma / (3.14 * D ** 2 * (1 - Arfg[ii-1]) * Phol / 4)
    cal_N2(Pc[ii-1])
    Vc[ii] = Ma / (3.14 * D ** 2 / 4 * Phol)
    Rel = Vc[ii] * D * Phol / Dvisl
    Prl = Dvisl * Cpl / Lbdl
    Reg = Vc[ii] * D * Phog / Dvisg
    Prg = Dvisg * Cpg / Lbdg            
    ReDB = 4 * Ma / (3.14 * D * Dvisl)

    if  Tw[ii] >= Ts:  
        if Y[ii-1] > 0.001 :
            xtt = ((1 - Y[ii-1]) / Y[ii-1]) ** 0.9 * (Phog / Phol) ** 0.5 * (Dvisl / Dvisg) ** 0.1
            F0 = (1 / xtt + 0.213) ** 0.736
        else:
            xtt = ((1 - 0.001) / 0.001) ** 0.9 * (Phog / Phol) ** 0.5 * (Dvisl / Dvisg) ** 0.1
            F0 = (1 / xtt + 0.213) ** 0.736
        Hfc = 0.023 * Rel ** 0.8 * Prl ** 0.4 * Lbdl / D
        ps1 = cal_Ps(Ts) * 100000 
        ps2 = cal_Ps(Tw[ii]) * 100000
        Hpb = 0.00122 * (Lbdl ** 0.79 * Cpl ** 0.45 * Phol ** 0.49 / (Srften ** 0.5 * Dvisl ** 0.29 * latent ** 0.24 * Phog ** 0.24)) * abs(Tw[ii] - Ts) ** 0.24 * abs(ps2 - ps1) ** 0.75
        REp = Rel * F0 ** 1.25 * 0.0001
        s = (1 + 0.00000253 * REp ** 1.17 * F0 ** 1.4625) ** (-1)
        Hchen = Hfc * F0 + Hpb * s
        h[ii] = Hchen

        hDB = 0.023 * ReDB ** 0.8 * Prl ** 0.4 * (Lbdl / D)
        if h[ii] < hDB: h[ii] = hDB; 
    else: 
        h[ii] = 0.023 * ReDB ** 0.8 * Prl ** 0.4 * (Lbdl / D)
    
    deta_T = (h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx) / (Ma * Cpl)
    deta_H = (h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx) / Ma

    if Rel > 2000 and Rel <= 50000 :
        ploss = 0.3164 / Rel ** 0.25 * dx * Phol * Vc[ii] ** 2 / (2 * D) 
    elif Rel < 2000 :
        ploss = 0.5 * 64 / Rel * Phol * Vc[ii] ** 2 * (dx / D)
    else:
        ploss = 0.184 / Rel ** 0.2 * dx * Phol * Vc[ii] ** 2 / (2 * D)
    Pc[ii] = (Pc[ii-1] * 100000 - ploss - (Vc[ii] ** 2 - V_temp ** 2) * Phol / 2) / 100000
    Tc[ii] = Tc[ii-1] + deta_T
    Hc[ii] = Hc[ii-1] + deta_H

    Y[ii] = Y[ii-1]
    Arfg[ii] = Arfg[ii-1]
    
    cal_N2(Pc[ii])
    if Tc[ii] > Ts:
        Tc[ii] = Ts
        dm = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx)) / latent
        if dm < 0 : dm = 0.000001
        Mg = Ma * Y[ii-1] + dm
        Ml = Ma - Mg
        Y[ii] = Mg / (Ml + Mg)
        Arfg[ii] = 1 - (Ml * Phog / Mg) / (Phol + Ml * Phog / Mg)
        Vc[ii] = Ml / (3.14 * D ** 2 * (1 - Arfg[ii]) * Phol / 4)
        if Arfg[ii] >= 1:
            Y[ii] = 1; Arfg[ii] = 1
            Reg = Mg * D / (Dvisg * 3.14 * D ** 2 / 4); Prg = Dvisg * Cpg / Lbdg; h[ii] = 0.023 * Reg ** 0.8 * Prg ** 0.4 * Lbdg / D
            deta_T = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx) - (Ml * Cpl * (Tc[ii] - Tc[ii-1])) - (Mg * Cpg * (Tc[ii] - Tc[ii-1]))) / (Ma * Cpg)
            deta_H = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx)) / Ma
            Tc[ii] = Tc[ii-1] + deta_T; Hc[ii] = Hc[ii-1] + deta_H; Vc[ii] = Ma / (3.14 * D ** 2 * Phog / 4)
            Mg = Ma; Ml = 0

def cal_phase(ii):
    global Ma, Ml, Mg, TONB, TFL, qCHF, boil_index
    Mg = Ma * Y[ii-1]
    Ml = Ma - Mg   
    Vc[ii] = Ml / (3.14 * D ** 2 * (1 - Arfg[ii-1]) * Phol / 4)
    Prg = Dvisg * Cpg / Lbdg
    Prl = Dvisl * Cpl / Lbdl
    Rel = Phol * Vc[ii] * D / Dvisl
    Reg = Phog * Vc[ii] * D / Dvisg
    Re = 4 * Ma / (3.14 * D * Dvisg)
    ReDB = 4 * Ma / (3.14 * D * Dvisl)
    xe = 0.6

    if Tw[ii] >= TFL : 
        Htc = 0.112 * (1 + (D / (Lc[ii] + 0.1 * D)) ** 0.7) * (Lbdg / D) * Re ** 0.57 * Prg ** 0.4 * (1 - xe) ** (-0.8)
        Hb = 0.05 * ((Phog * (Phol - Phog) * 9.8 * Lbdg ** 3 * (latent + 0.375 * Cpg * (Tw[ii] - Ts))) / (4 * Lc[ii] * Dvisg * (Tw[ii] - Ts))) ** 0.25
        Hdw = 0.00007 * (Lbdl / D) * math.exp(-((Lc[ii] / D) ** 0.45)) * ((Cpl * (Tw[ii] - Ts)) / latent) ** (-1) * ((1.62 * Ml ** 2) / (D ** 3 * Phol * Srften)) ** 2
        h[ii] = Htc + Hb + Hdw
        boil_index = 1
    elif Tw[ii] <= Ts :   
        h[ii] = 0.023 * ReDB ** 0.8 * Prl ** 0.4 * (Lbdl / D)
    else: 
        xtt = ((1 - Y[ii-1]) / Y[ii-1]) ** 0.9 * (Phog / Phol) ** 0.5 * (Dvisl / Dvisg) ** 0.1
        F0 = (1 / xtt + 0.213) ** 0.736
        Hfc = 0.023 * Rel ** 0.8 * Prl ** 0.4 * Lbdl / D
        ps1 = cal_Ps(Ts) * 100000
        ps2 = cal_Ps(Tw[ii]) * 100000
        Hpb = 0.00122 * (Lbdl ** 0.79 * Cpl ** 0.45 * Phol ** 0.49 / (Srften ** 0.5 * Dvisl ** 0.29 * latent ** 0.24 * Phog ** 0.24)) * abs(Tw[ii] - Ts) ** 0.24 * abs(ps2 - ps1) ** 0.75
        REp = Rel * F0 ** 1.25 * 0.0001
        s = (1 + 0.00000253 * REp ** 1.17 * F0 ** 1.4625) ** (-1)
        Hchen = Hfc * F0 + Hpb * s
        qchen = Hchen * (Tw[ii] - Tc[ii-1])
    
    if boil_index == 1:
        h[ii] = 83
        
    if Rel < 2000: n = 1;cl = 64
    elif Rel >= 2000 and Rel < 5000: n = 0.25;cl = 0.316
    else: n = 0.2;cl = 0.184
    if Reg < 2000: m = 1;cl = 64
    elif Reg >= 2000 and Reg < 50000: m = 0.25;cl = 0.316
    else: m = 0.2;cl = 0.184

    if Rel < 2000 and Reg < 2000 : C = 5
    elif Rel >= 2000 and Reg < 2000 : C = 10
    elif Rel < 2000 and Reg > 2000 : C = 12
    else: C = 20

    Y[ii] = Mg / (Ml + Mg)
    X = (Phog / Phol * (1 - Y[ii]) ** 1.8 / Y[ii] ** 1.8 * (Dvisl / Dvisg) ** 0.2) ** 0.5
    phil = (X ** 2 + C * X + 1) ** 0.5 / X
    if Rel < 5000 and Reg > 50000 : n = 0.2; cl = 0.184; phil = (Phol / Phog) ** 0.5
    F = cl * Rel ** (-n)
    G = Ma / (3.14 * D ** 2 / 4)
    ploss = (2 * F * dx * G ** 2 * (1 - Y[ii]) ** 2 / (4 * D * Phol)) * phil ** 2 
    deta_H = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx)) / Ma

    cal_N2(Pc[ii-1])
    Tc[ii] = Ts
    Hc[ii] = Hc[ii-1] + deta_H
    dm = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx)) / latent
    if dm < 0 : dm = 0.0
    Mg = Mg + dm
    Ml = Ma - Mg
    Y[ii] = Mg / (Ml + Mg)
    Arfg[ii] = 1 - (Ml * Phog / Mg) / (Phol + Ml * Phog / Mg)
    Vc[ii] = Ml / (3.14 * D ** 2 * (1 - Arfg[ii]) * Phol / 4)
    Pc[ii] = (Pc[ii-1] * 100000 - ploss - Vc[ii] ** 2 * (Phol * (1 - Arfg[ii]) + Arfg[ii] * Phog) / 2 + Vc[ii] ** 2 * (Phol * (1 - Arfg[ii-1]) + Arfg[ii-1] * Phog) / 2) / 100000
    if Arfg[ii] >= 1 :
        Y[ii] = 1;Arfg[ii] = 1
        Reg = Mg * D / (Dvisg * 3.14 * D ** 2 / 4);Prg = Dvisg * Cpg / Lbdg; h[ii] = 0.023 * Reg ** 0.8 * Prg ** 0.4 * Lbdg / D
        deta_T = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx) - (Ml * Cpl * (Tc[ii] - Tc[ii-1])) - (Mg * Cpg * (Tc[ii] - Tc[ii-1]))) / (Ma * Cpg)
        deta_H = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx)) / Ma
        Tc[ii] = Tc[ii-1] + deta_T; Hc[ii] = Hc[ii-1] + deta_H; Vc[ii] = Ma / (3.14 * D ** 2 * Phog / 4)
        Mg = Ma; Ml = 0

def cal_g(ii):
    if Tc[ii-1] < 120 :
        Cpg = (2.891738674 - 0.057984266 * Tc[ii-1] + 0.000707606 * Tc[ii-1] ** 2 - 0.00000392931 * Tc[ii-1] ** 3 + 0.00000000831 * Tc[ii-1] ** 4) * 1000
        Lbdg = (0.106985591 + 0.074993969 * Tc[ii-1] + 0.000376941 * Tc[ii-1] ** 2 - 0.00000253364 * Tc[ii-1] ** 3 + 0.00000000538422 * Tc[ii-1] ** 4) * 0.001
        Dvisg = (0.026978089 + 0.068875284 * Tc[ii-1] + 0.0000630018 * Tc[ii-1] ** 2 - 0.000000755697 * Tc[ii-1] ** 3 + 0.00000000169223 * Tc[ii-1] ** 4) * 0.000001
        Phog = 20.67603989 - 0.450301075 * Tc[ii-1] + 0.00466498 * Tc[ii-1] ** 2 - 0.0000235267 * Tc[ii-1] ** 3 + 0.0000000466757 * Tc[ii-1] ** 4
    elif Tc[ii-1] <= 300 :
        Phog = (Pc[ii-1] * 100000 * 28) / (8314 * Tc[ii-1])
        Cpg = (1.209680332 - 0.002660496 * Tc[ii-1] + 0.0000162696 * Tc[ii-1] ** 2 - 0.0000000450562 * Tc[ii-1] ** 3 + 4.71909E-11 * Tc[ii-1] ** 4) * 1000
        Dvisg = (-0.417307254 + 0.082584193 * Tc[ii-1] - 0.000098396 * Tc[ii-1] ** 2 + 0.000000103259 * Tc[ii-1] ** 3 - 4.94201E-11 * Tc[ii-1] ** 4) * 0.000001
        Lbdg = (-1.155836746 + 0.114783343 * Tc[ii-1] - 0.00010168 * Tc[ii-1] ** 2 + 0.0000000683606 * Tc[ii-1] ** 3 - 6.35283E-13 * Tc[ii-1] ** 4) * 0.001
    else:
        Phog = (Pc[ii-1] * 100000 * 28) / (8314 * Tc[ii-1])
        Cpg = (1.027977283 + 0.000242958 * Tc[ii-1] - 0.00000132728 * Tc[ii-1] ** 2 + 0.00000000268404 * Tc[ii-1] ** 3 - 1.54763E-12 * Tc[ii-1] ** 4) * 1000
        Dvisg = (0.107894047 + 0.076774811 * Tc[ii-1] - 0.0000746043 * Tc[ii-1] ** 2 + 0.0000000609547 * Tc[ii-1] ** 3 - 2.24127E-11 * Tc[ii-1] ** 4) * 0.000001
        Lbdg = (-0.778133404 + 0.11169081 * Tc[ii-1] - 0.0000958852 * Tc[ii-1] ** 2 + 0.0000000777002 * Tc[ii-1] ** 3 - 0.000000000028246 * Tc[ii-1] ** 4) * 0.001

    Reg = Ma * D / (Dvisg * 3.14 * D ** 2 / 4)
    Prg = Dvisg * Cpg / Lbdg
    if Tw[ii] < (Tc[ii-1] + 50) :
        h[ii] = 0.023 * Reg ** 0.8 * Prg ** 0.4 * Lbdg / D
    elif Tw[ii] >= (Tc[ii-1] + 50) :
        h[ii] = 0.023 * Reg ** 0.8 * Prg ** 0.4 * (Tc[ii-1] / Tw[ii]) ** 0.5 * Lbdg / D

    deta_T = (h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx) / (Ma * Cpg)
    deta_H = (h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx) / Ma
    Tc[ii] = Tc[ii-1] + deta_T
    Hc[ii] = Hc[ii-1] + deta_H
    Y[ii] = Y[ii-1]
    Arfg[ii] = Arfg[ii-1]    
    Phog = (Pc[ii-1] * 100000 * 28) / (8314 * Tc[ii])
    Vc[ii] = Ma / (3.14 * D ** 2 / 4 * Phog)
    if Reg > 2000 and Reg <= 50000 :
        ploss = 0.3164 / Reg ** 0.25 * dx * Phog * Vc[ii] ** 2 / (2 * D)
    elif Reg < 2000 :
        ploss = 64 / Reg * dx * Phog * Vc[ii] ** 2 / (2 * D)
    else:
        ploss = 0.184 / Reg ** 0.2 * dx * Phog * Vc[ii] ** 2 / (2 * D)
    Pc[ii] = (Pc[ii-1] * 100000 - ploss - (Vc[ii] ** 2 - Vc[ii-1] ** 2) * Phog / 2) / 100000
    
    cal_N2(Pc[ii])
    if Tc[ii] < Ts:
        Tc[ii] = Ts
        dm = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx)) / latent
        Mg = Ma * Y[ii-1] + dm
        Ml = Ma - Mg
        Y[ii] = Mg / (Ml + Mg)
        Arfg[ii] = 1 - (Ml * Phog / Mg) / (Phol + Ml * Phog / Mg)
        Vc[ii] = Ml / (3.14 * D ** 2 * (1 - Arfg[ii]) * Phol / 4)
        if Arfg[ii] >= 1:
            Y[ii] = 1; Arfg[ii] = 1
            Reg = Mg * D / (Dvisg * 3.14 * D ** 2 / 4); Prg = Dvisg * Cpg / Lbdg; h[ii] = 0.023 * Reg ** 0.8 * Prg ** 0.4 * Lbdg / D
            deta_T = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx) - (Ml * Cpl * (Tc[ii] - Tc[ii-1])) - (Mg * Cpg * (Tc[ii] - Tc[ii-1]))) / (Ma * Cpg)
            deta_H = ((h[ii] * (Tw[ii] - Tc[ii-1])- hlng * (Tc[ii-1] - Tinlet)) * (3.14 * D * dx)) / Ma
            Tc[ii] = Tc[ii-1] + deta_T; Hc[ii] = Hc[ii-1] + deta_H; Vc[ii] = Ma / (3.14 * D ** 2 * Phog / 4)
            Mg = Ma; Ml = 0

def cal_T_profile():
    global n, dx, dt, Dcu
    a = np.zeros(n+1); b = np.zeros(n+1); c = np.zeros(n+1); d = np.zeros(n+1)
    for i in range(n+1):
        cal_Cu(Tw[i])
        b[i] = - Lbdcu *dt / dx; c[i] = - Lbdcu *dt / dx
        d[i] = 4 / Dcu * (h[i]+hLNG) * Tc[i] * dx * dt + Phocu * Cpcu * dx * Tw[i] + Source[i] * dx * dt + Qp[i]*4/Dcu * dx * dt
        a[i] = 4 / Dcu * (h[i]+hLNG) * dx * dt + Phocu * Cpcu * dx + 2 * Lbdcu * dt / dx
        if i == 0 or i == n :
            a[i] = a[i] / 2; d[i] = d[i] / 2
    c[0] = 0; b[n] = 0
    p = [0]*(n+1); q = [0]*(n+1)
    for i in range(n+1):
        if not i:
            p[i] = b[i] / a[i]
            q[i] = d[i] / a[i]
        else:
            p[i] = b[i] / (a[i] - p[i-1] * c[i])
            q[i] = (d[i] - q[i-1] * c[i]) / (a[i] - p[i-1] * c[i])
    for i in range(n, -1, -1):
        if i == n:
            Tw[i] = q[i]
        else:
            Tw[i] = q[i] - p[i] * Tw[i+1]

def cal_current(ii):
    cal_Cu(Tw[ii])
    miu0 = 4*math.pi*0.0000001; Pi = 0.15; aiaj = 1; index = 32; E0 = 0.0001
    Ltap = miu0 * math.pi * (D/2/Pi) ** 2 + miu0 * math.log(D/Dcu) /(2*math.pi)
    Mtap = math.pi * (Dcu/2)**2 * aiaj * miu0 / Pi**2 + miu0 * math.log(D/Dcu) / (2*math.pi)
    Lcu = miu0 / (8*math.pi) + miu0 * math.log(D/Dcu) /(2*math.pi)
    Mcu = miu0 * math.log(D/Dcu) /(2*math.pi)
    Ic[ii] = 690 - 7.5 * Tw[ii]
    if Ic[ii] > 0:
        a=np.zeros(index+1)
        a[0]=E0/Ic[ii]**index
        a[index-1]=m*Rcu*4/(math.pi*Dcu**2) + Ltap/dt - (m-1)*Mtap/dt + 2*m*Mcu/dt + m*Lcu/dt
        a[index]=-Itap[ii]*(Ltap/dt - (m-1)*Mtap/dt + 2*m*Mcu/dt + m*Lcu/dt) - Rcu*4/(math.pi*Dcu**2)*Iop
        result=np.roots(a)
        for j in range(len(result)):
            if np.imag(result[j]) == 0: Itap[ii] = np.real(result[j])
    else:
        Itap[ii] = 0
    Icu[ii] = Iop - m*Itap[ii]
    Source[ii] = Rcu * Icu[ii]**2 / (math.pi*Dcu**2/4)**2

class ProgressBar():
    def __init__(self, width=50):
        self.width = width
    def update(self, x, t):
        assert 0 <= x <= 100 
        pointer = int(self.width * (x / 100.0))
        sys.stdout.write( '\rCalculate Time%.2fs/%ds [%s]%d%% Total Time%.0fs' %(t, totalT, '#' * pointer + '.' * (self.width - pointer), int(x), count2-count1))
        sys.stdout.flush()
        if x == 100: print ('')

file_path = os.getcwd()
f1 = []; record_time = [1,2,5,10,15,20,25,30,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]
for i in range(len(record_time)):
    f1.append(open(r"%s\t=%s.txt"%(file_path, record_time[i]), "w"))
f2 = []; record_point = [0.40,0.45,0.50,0.55,0.60]
for i in range(len(record_point)):
    f2.append(open(r"%s\x=%.2f.txt"%(file_path, record_point[i]), "w"))

Tw = np.zeros(n+1); h = np.zeros(n+1); boil = list(range(n+1))
Pc = np.zeros(n+1); Tc = np.zeros(n+1); Vc = np.zeros(n+1); Hc = np.zeros(n+1); Arfg = np.zeros(n+1); Y = np.zeros(n+1)
Lc = np.zeros(n+1); Itap = np.zeros(n+1); Icu = np.zeros(n+1); Ic = np.zeros(n+1); Source = np.zeros(n+1)
for i in range(n+1):
    Tw[i] = 80
    Itap[i] = Iop / m
    if i >= 450 and i <= 550:
        Qp[i] = 1000/(1*math.pi*Dcu*0.1)  
t = 0; Qwall = 0
cal_TF()
    
for i in range(len(record_time)): f1[i].close()
for i in range(len(record_point)):f2[i].close()
