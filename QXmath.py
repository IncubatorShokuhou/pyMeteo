import math


def showalter(t8, td8, t5):  # 输入850hPa温度和露点温度，500hPa温度，计算沙氏指数。单位：摄氏度
    P5 = 500
    P8 = 850
    Cpd = 0.2403
    Cpv = 0.445
    Rd = 6.85578*0.01
    Rw = 11.017874*0.01
    T0 = 273.15
    L0 = 597.40
    C = 1.002
    C1 = 0.57
    step = 10
    T8 = T0+t8
    Etd8 = E_WATER(t8)
    w = (Rd/Rw)*Etd8/(P8-Etd8)
    ml = (Cpd*(1+Cpv*w/Cpd))/(Rd*(1+w/(Rd/Rw)))
    Ta = Tc(P8, t8, td8)+T0
    Pa = P8*(pow(Ta/T8, ml))
    m2 = (Cpd/Rd)*(1+C*w/Cpd)
    ETa = E_WATER(Ta-T0)
    Qc = math.log((Pa-ETa)/pow(Ta, m2))-(0.622/Rd) * \
        (((L0+C1*(T0-Ta))/Ta)*(ETa/(Pa-ETa)))
    out = T8
    Eout = E_WATER(out-T0)
    Q5 = math.log((P5-Eout)/pow(out, m2))-(0.622/Rd) * \
        (((L0+C1*(T0-out))/out)*(Eout/(P5-Eout)))
    while (abs(Qc-Q5) > 0.0001):
        if (Qc > Q5):
            out -= step
        else:
            out += step
            step /= 5
            out -= step
        Eout = E_WATER(out-T0)
        Q5 = math.log((P5-Eout)/pow(out, m2))-(0.622/Rd) * \
            (((L0+C1*(T0-out))/out)*(Eout/(P5-Eout)))
    return t5-(out-T0)


def E_WATER(td):  # 计算水面饱和蒸汽压
    E0 = 6.1078
    T0 = 273.15
    Cl = 0.57
    Rw = 0.1101787372
    L0 = 597.4
    T = T0+td
    E = ((L0+Cl*T0)*(T-T0))/(Rw*T0*T)
    E = math.exp(E)
    E = E*E0*pow(T0/T, Cl/Rw)
    return E


def Tc(p, t, td):  # 根据温度和露点计算凝结高度上的温度
    step = 10.0
    Cpd = 0.2403
    Cpv = 0.445
    Rd = 6.85578*0.01
    Rw = 11.017874*0.01
    T0 = 273.15
    T = T0+t
    Td = T0+td
    Etd = E_WATER(td)
    w = (Rd/Rw)*Etd/(p-Etd)
    ml = (Cpd*(1+Cpv*w/Cpd))/(Rd*(1+w/(Rd/Rw)))
    Z0 = pow(T, ml)/Etd
    out = Td
    Z = pow(out, ml)/E_WATER(out-T0)
    while abs(Z-Z0) > 10:
        if (Z < Z0):
            out = out-step
        else:
            out = out+step
            step = step/5
            out = out-step
        Z = pow(out, ml)/E_WATER(out-T0)
        if (step < 0.00000001):
            break
    return out-T0


def K(T850, Td850, T700, Td700, T500):  # K指数
    return T850-T500+Td850-(T700-Td700)


def A(T850, Td850, T700, Td700, T500, Td500):  # A指数
    return T850-T500-(T850-Td850)-(T700-Td700)-(T500-Td500)


def ttd700(T700, Td700):  # 700hPa温度露点差
    return T700-Td700


def ttd850(T850, Td850):  # 850hPa温度露点差
    return T850-Td850


def ttd925(T925, Td925):  # 925hPa温度露点差
    return T925-Td925


def tt500(T850, T500):  # 大概是温度差？？
    return T850-T500


def thunder_potential(T925, Td925, T850, Td850, T700, Td700, T500, Td500):
    if (K(T850-273.15, Td850-273.15, T700-273.15, Td700-273.15, T500-273.15) > 33.0):
        k1 = 1
    else:
        k1 = 0

    if (showalter(T850-273.15, Td850-273.15, T500-273.15) < 0.0):
        k2 = 1
    else:
        k2 = 0

    if (A(T850-273.15, Td850-273.15, T700-273.15, Td700-273.15, T500-273.15, Td500-273.15) > 10.0):
        k3 = 1
    else:
        k3 = 0

    if (ttd700(T700, Td700) < 3.0):
        k4 = 1
    else:
        k4 = 0

    if (ttd850(T850, Td850) < 3.0):
        k5 = 1
    else:
        k5 = 0

    if (ttd925(T925, Td925) < 3.0):
        k6 = 1
    else:
        k6 = 0

    if (tt500(T850, T500) > 23.0):
        k7 = 1
    else:
        k7 = 0

    y = 0.110+0.154*k1+0.142*k2+0.061*k3+0.146*k4+0.131*k5+0.03*k6+0.097*k7

    if (y > 0.5):
        k8 = 1
    else:
        k8 = 0
    return k8

def relhum_ttd(T,Td):   # 根据温度和露点计算相对湿度
    gc  = 461.5             # [j/{kg-k}]   gas constant water vapor
    gc  = gc/(1000.*4.186)  # [cal/{g-k}]  change units
                                       # lhv=latent heat vap
    lhv = ( 597.3-0.57*(T-273.15) )        # dutton top p273 [empirical]

    rh  = math.exp( (lhv/gc)*(1.0/T - 1.0/Td) ) 
    #rh=math.exp(lhv/gc/T)/math.exp(lhv/gc/Td)
    rh  = rh*100.0
    return rh

def dewtemp_trh(Tk,RH):   # 根据温度和相对湿度计算露点
    '''
    calculate the dew pt temperature [k]:
    input:
        rh - relative humidty [%]
        tk - temperature [k]
    output:
        tdk- dew point temperature
             equation used is Dutton p274 (6) [k]

    local:
        gc  - gas constant for water vapor [j/{kg-k}]
              gcx [cal/{g-k}]
        lhv - latent heat vap (Dutton top p273)
    '''
    GC=461.5
    GCX=GC/ (1000*4.186)
    LHV = (597.3-0.57* (Tk-273.15))/GCX
    TDK = Tk*LHV/ (LHV-Tk*math.log(RH*0.01))
    return TDK

if __name__ == "__main__":
    #a = showalter(16.6, 0.6, -15.9)
    tk =   18. + 273.15   # K
    rh =   46.5           # %

    td = dewtemp_trh(tk,rh) - 273.15  # td = 6.3 C
    print(td)