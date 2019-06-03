import math
def showalter(t8,td8,t5):  #输入850hPa温度和露点温度，500hPa温度，计算沙氏指数。单位：摄氏度
    P5=500
    P8=850
    Cpd=0.2403
    Cpv=0.445
    Rd=6.85578*0.01
    Rw=11.017874*0.01
    T0=273.15
    L0=597.40
    C=1.002
    C1=0.57
    step=10
    T8=T0+t8
    Etd8=E_WATER(t8)
    w=(Rd/Rw)*Etd8/(P8-Etd8)
    ml=(Cpd*(1+Cpv*w/Cpd))/(Rd*(1+w/(Rd/Rw)))
    Ta=Tc(P8,t8,td8)+T0
    Pa=P8*(pow(Ta/T8,ml))
    m2=(Cpd/Rd)*(1+C*w/Cpd)
    ETa=E_WATER(Ta-T0)
    Qc=math.log((Pa-ETa)/pow(Ta,m2))-(0.622/Rd)*(((L0+C1*(T0-Ta))/Ta)*(ETa/(Pa-ETa)))
    out=T8
    Eout=E_WATER(out-T0)
    Q5=math.log((P5-Eout)/pow(out,m2))-(0.622/Rd)*(((L0+C1*(T0-out))/out)*(Eout/(P5-Eout)))
    while (abs(Qc-Q5)>0.0001):
        if(Qc>Q5):
            out-=step
        else:
            out+=step
            step/=5
            out-=step
        Eout=E_WATER(out-T0)
        Q5=math.log((P5-Eout)/pow(out,m2))-(0.622/Rd)*(((L0+C1*(T0-out))/out)*(Eout/(P5-Eout)))
    return t5-(out-T0)

def E_WATER(td):   #计算水面饱和蒸汽压
    E0=6.1078
    T0=273.15
    Cl=0.57
    Rw=0.1101787372
    L0=597.4
    T=T0+td
    E=((L0+Cl*T0)*(T-T0))/(Rw*T0*T)
    E=math.exp(E)
    E=E*E0*pow(T0/T,Cl/Rw)
    return E

def Tc(p,t,td):   #根据温度和露点计算凝结高度上的温度
    step=10.0
    Cpd=0.2403
    Cpv=0.445
    Rd=6.85578*0.01
    Rw=11.017874*0.01
    T0=273.15
    T=T0+t
    Td=T0+td
    Etd=E_WATER(td)
    w=(Rd/Rw)*Etd/(p-Etd)
    ml=(Cpd*(1+Cpv*w/Cpd))/(Rd*(1+w/(Rd/Rw)))
    Z0=pow(T,ml)/Etd
    out=Td
    Z=pow(out,ml)/E_WATER(out-T0)
    while abs(Z-Z0) > 10:
        if (Z<Z0):
            out=out-step
        else:
            out=out+step
            step=step/5
            out=out-step
        Z=pow(out,ml)/E_WATER(out-T0)
        if (step < 0.00000001): break
    return out-T0


if __name__ == "__main__":
    a=showalter(16.6,0.6,-15.9)
    print(a)