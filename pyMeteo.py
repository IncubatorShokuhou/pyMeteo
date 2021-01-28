import math
import numpy as np
from geopy import distance

@np.vectorize
def showalter(t8, td8, t5):
    # 输入850hPa温度和露点温度，500hPa温度，计算沙氏指数。单位：摄氏度
    '''
    caculating Showalter Index 
    Input:
        t8: temperature at 850hPa [degree cetigrade ]
        td: dew point temperature at 850hPa [degree cetigrade ] 
        t5: temperature at 850hPa [degree cetigrade ]
    Output:
        showalter index
    '''

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

@np.vectorize
def E_WATER(td):  
    '''
    caculate aturated vapor pressure
    Input:
        td: dew point temperature [degree cetigrade] 
    Output:
        showalter index

    '''
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

@np.vectorize
def Tc(p, t, td):
    '''
    Calculate the temperature at condensation height based on temperature and dew point
    Input:
        p: pressure [hPa]
        t: temperature [degree cetigrade]
        td: dew point temperature [degree cetigrade]
    Output:
        temperature at condensation height [degree cetigrade]
    '''
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

@np.vectorize
def K(T850, Td850, T700, Td700, T500):  # K指数
    '''
    caculating K index
    The larger the K index, the more unstable the stratification. 
    The statistical results: 
        K<20 no thunderstorm; 
        20<K<25 isolated thunderstorm; 
        25<K<30 sporadic thunderstorm; 
        30<K<35 flakes of thunderstorm.
    Input:
        T850: temperature at 850hPa [degree cetigrade]
        Td850: dew point temperature at 850hPa [degree cetigrade]
        T700: temperature at 700hPa [degree cetigrade]
        Td700: dew point temperature at 700hPa [degree cetigrade]
        T500: temperature at 500hPa [degree cetigrade]
    Output:
        K index
    '''
    return T850-T500+Td850-(T700-Td700)

@np.vectorize
def A(T850, Td850, T700, Td700, T500, Td500):  
    '''
    caculating A index
    Input:
        T850: temperature at 850hPa [degree cetigrade]
        Td850: dew point temperature at 850hPa [degree cetigrade]
        T700: temperature at 700hPa [degree cetigrade]
        Td700: dew point temperature at 700hPa [degree cetigrade]
        T500: temperature at 500hPa [degree cetigrade]
        Td500: dew point temperature at 500hPa [degree cetigrade]
    Output:
        A index
    '''
    return T850-T500-(T850-Td850)-(T700-Td700)-(T500-Td500)

@np.vectorize
def ttd700(T700, Td700):  
    # Temperature-dew point difference at 700hPa
    return T700-Td700

@np.vectorize
def ttd850(T850, Td850):  
    # Temperature-dew point difference at 850hPa
    return T850-Td850

@np.vectorize
def ttd925(T925, Td925):  
    # Temperature-dew point difference dat 925hhPa
    return T925-Td925

@np.vectorize
def tt500(T850, T500):  
    # Temperature difference bewteen 850hPa and 500hPa
    return T850-T500

@np.vectorize
def relhum_ttd(T,Td):   
    '''
    Calculate relative humidity giving temperature and dew point
    Input:
        T: temperature [K]
        Td: dew point temperature [K]
    Output:
        relative humidity [%]
    '''
    gc  = 461.5             # [j/{kg-k}]   gas constant water vapor
    gc  = gc/(1000.*4.186)  # [cal/{g-k}]  change units
                                       # lhv=latent heat vap
    lhv = ( 597.3-0.57*(T-273.15) )        # dutton top p273 [empirical]
    rh  = math.exp( (lhv/gc)*(1.0/T - 1.0/Td) ) 
    #rh=math.exp(lhv/gc/T)/math.exp(lhv/gc/Td)
    rh  = rh*100.0
    return rh

@np.vectorize
def dewtemp_trh(Tk,RH):  
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

@np.vectorize
def TT(T850,R850,T500):  
    '''
    total totals index
    Input: 
        T850: temperature at 850hPa [K]
        R850: relative humidty at 850hPa [%]
        T500: temperature at 500hPa [K]
    '''
    Td850 = dewtemp_trh(T850,R850)
    TT = T850 + Td850 - 2*T500
    return TT

@np.vectorize
def ws(u,v): 
    '''
    caculating wind speed giving u and v
    Input:
        u: u-wind
        v: v-wind
    Output:
        wind speed
    '''
    return math.sqrt(u**2+v**2)

@np.vectorize
def wd(u,v): 
    '''
    caculating wind direction giving u and v
    Input:
        u: u-wind
        v: v-wind
    Output:
        wind direction [degree]
    '''
    return 180+math.atan2(u,v)*180/math.pi

@np.vectorize
def u(ws,wd): 
    '''
    caculating u-wind giving wind speed and wind direction
    Input:
        ws: wind speed
        wd: wind direction [degree]
    Output:
        u-wind
    '''
    return -ws*math.sin(wd/180*math.pi)
  
@np.vectorize
def v(ws,wd): 
    '''
    caculating v-wind giving wind speed and wind direction
    Input:
        ws: wind speed
        wd: wind direction [degree]
    Output:
        v-wind
    '''
    return -ws*math.cos(wd/180*math.pi)
    
@np.vectorize
def SWEAT_caculate(T850,R850,T500,U850,V850,U500,V500):
    '''
    caculating SWEAT Severe Weather Threat Index
    SWEAT=12*Td850 + 20*(TT-49) + 4*WS850 + 2*WS500 + 125*(sin(WD500-WD850)+0.2)
    Input:
        T850: tempereature at 850hPa [K]
        R850: relative humidty at 850hPa [%]
        T500: tempereature at 500hPa [K]
        U850: u-wind at 850hPa [m/s]
        V850: v-wind at 850hPa [m/s]
        U500: u-wind at 500hPa [m/s]
        V500: v-wind at 500hPa [m/s]
    Output:
        SWEAT Index
    '''
    S_1=12*dewtemp_trh(T850,R850)
    S_2=20*(TT(T850,R850,T500)-49)
    S_3=4*ws(U850,V850)
    S_4=2*ws(U500,V500)
    S_5=125*(math.sin((wd(U500,V500)-wd(U850,V850))/180*math.pi)+0.2)
    if S_1 <0 : S_1 = 0
    if S_2 <0 : S_2 = 0
    if S_3 <0 : S_3 = 0
    if S_4 <0 : S_4 = 0
    if S_5 <0 : S_5 = 0
    if wd(U850,V850) < 130 or wd(U850,V850) > 250 or wd(U500,V500) < 210 or wd(U500,V500) > 310 : S_5 = 0
    return S_1+S_2+S_3+S_4+S_5

@np.vectorize
def earth_distance(lat1, lon1, lat2, lon2,unit="km"):
    '''
        caculating the distance bewettn two point on the earth giving the latitude and longiude of the two points
        Input:
            lat1: latitude of point 1
            lon1: longiude of point 1
            lat2: latitude of point 2
            lon2: longiude of point 2
            unit: unit of distance. can be "m","meter","km" ...
        Output:
            distance, with certain unit giving in the parameter
    '''
    a=distance.distance((lat1, lon1), (lat2, lon2))
    return eval("a."+unit)

@np.vectorize
def relhum(T,W,P):
    '''
    caculating relative humidity using temperature,mixing ratio and pressure
    Input:
        T:  temperature  [K]
        W:  mixing ratio [kg/kg]
        P:  pressure     [Pa]
    Output:
        relative humidity [%]
    '''
    TABLE = [0.01403,0.01719,0.02101,0.02561,0.03117,0.03784,
            0.04584,0.05542,0.06685,0.08049,0.09672,0.1160,0.1388,
            0.1658,0.1977,0.2353,0.2796,0.3316,0.3925,0.4638,
            0.5472,0.6444,0.7577,0.8894,1.042,1.220,1.425,
            1.662,1.936,2.252,2.615,3.032,3.511,4.060,
            4.688,5.406,6.225,7.159,8.223,9.432,10.80,
            12.36,14.13,16.12,18.38,20.92,23.80,27.03,
            30.67,34.76,39.35,44.49,50.26,56.71,63.93,
            71.98,80.97,90.98,102.1,114.5,128.3,143.6,
            160.6,179.4,200.2,223.3,248.8,276.9,307.9,
            342.1,379.8,421.3,466.9,517.0,572.0,632.3,
            698.5,770.9,850.2,937.0,1032.0,1146.6,1272.0,1408.1,1556.7,1716.9,1890.3,2077.6,
         2279.6,2496.7,2729.8,2980.0,3247.8,3534.1,
         3839.8,4164.8,4510.5,4876.9,5265.1,5675.2,
         6107.8,6566.2,7054.7,7575.3,8129.4,8719.2,
         9346.5,10013.0,10722.0,11474.0,12272.0,13119.0,
         14017.0,14969.0,15977.0,17044.0,18173.0,19367.0,
         20630.0,21964.0,23373.0,24861.0,26430.0,28086.0,
         29831.0,31671.0,33608.0,35649.0,37796.0,40055.0,
         42430.0,44927.0,47551.0,50307.0,53200.0,56236.0,
         59422.0,62762.0,66264.0,69934.0,73777.0,77802.0,
         82015.0,86423.0,91034.0,95855.0,100890.0,106160.0,
         111660.0,117400.0,123400.0,129650.0,136170.0,142980.0,
         150070.0,157460.0,165160.0,173180.0,181530.0,190220.0,
         199260.0,208670.0,218450.0,228610.0,239180.0,250160.0,
         261560.0,273400.0,285700.0,298450.0,311690.0,325420.0,
         339650.0,354410.0,369710.0,385560.0,401980.0,418980.0,
         436590.0,454810.0,473670.0,493170.0,513350.0,534220.0,
         555800.0,578090.0,601130.0,624940.0,649530.0,674920.0,
         701130.0,728190.0,756110.0,784920.0,814630.0,845280.0,
         876880.0,909450.0,943020.0,977610.0,1013250.0,
    #*PL*ERROR* Too many continuation lines generated
         1049940.0,1087740.0,1087740.0]
    
    TP = T
    if (TP > 375.16):
        TP = 375.16
    if (TP < 173.16):
        TP = 173.16
    IT = int(TP - 173.16)
    T2 = 173.16 + IT
    ES = (T2+1.-TP)*TABLE[IT] + (TP-T2)*TABLE[IT+1]
    ES = ES*0.1
    DRELHUM = (W* (P-0.378*ES)/ (0.622*ES))*100.0
    '''
    C There is some discussion about whether the humidities
    C should be allowed to be above 100%. Right now, we're
    C leaving the code as is. We may eventually comment
    C this line out, or create another function.
    C
    C Dennis decided to comment this line out for V5.2.0.
        C
    '''
    #   if (DRELHUM.GT.100.) DRELHUM = 100.
    if (DRELHUM < 0.0):
        DRELHUM = 0.0001
    return DRELHUM


@np.vectorize
def visibility(RH,T,method:str):
    '''
    estimating visibility using relativate humidity and temperature.
    Input:
        T:  temperature  [K]
        RH: relative humidty [%]
    Output:
        visibility [km]
    '''
    if method == "RUC" : # using Rapid Updata Cyclef method
        q = min(80.0,(RH/100-0.15))
        return 60.0 * math.exp(-2.5*q)  * 1000
    elif method == "FSL": # using Forecast System Laboratory method
        vis = 6000*(T-dewtemp_trh(T,RH))/(RH**1.75)
        return vis
def g(lat):
    '''
    caculate more accurate gravitational acceleration using latitude.
    '''
    return 9.7803*(1+0.0053024*(math.sin(lat))**2-0.000005*(math.sin(2*lat))**2)
if __name__ == "__main__":
    #a = showalter(16.6, 0.6, -15.9)
    tk =   18. + 273.15   # K
    rh =   46.5           # %

    td = dewtemp_trh(tk,rh) - 273.15  # td = 6.3 C
    print(td)
