# 说明
[English version](https://github.com/IncubatorShokuhou/pyMeteo/blob/master/README.md)   
===========================   
一些简单的气象诊断代码。参考了ncl源代码，`李社宏.用C语言开发的气象常用参数和物理量计算函数库（一）[J].陕西气象,1994(03):42-45.`，以及其他一些资料。
===========================   

一些函数的简单说明:

`showalter(t8,td8,t5)`:  #输入850hPa温度和露点温度，500hPa温度，计算沙氏指数。单位：摄氏度  PS:算沙氏指数的时候出bug,八成是R有负值。基本上都是网格化数据插值插出来的问题。

`E_WATER(td)`:   #计算水面饱和蒸汽压

`Tc(p,t,td)`:   #根据温度和露点计算凝结高度上的温度

`relhum_ttd(T,Td)`:   # 根据温度和露点计算相对湿度(T:temperature[K],Td:temperature[K])

`dewtemp_trh(Tk,RH)`:   # 根据温度和相对湿度计算露点

`K(T850, Td850, T700, Td700, T500)`:  # K指数

`A(T850, Td850, T700, Td700, T500, Td500)`:  # A指数

`earth_distance(lat1, lon1, lat2, lon2)`:  #两点间距离。单位KM

`TT(T850,R850,T500)`:  #TT全总指数

`ws(u,v)`: #风速

`wd(u,v)`: #风向，返回角度制

`u(ws,wd)`: #计算u风。输入风速和风向（角度制）

`v(ws,wd)`: #计算v风。输入风速和风向（角度制）

`SWEAT_calculate(T850,R850,T500,U850,V850,U500,V500)`:   # SWEAT天气强威胁指数

`relhum(T,W,P)`:  # 计算相对湿度(T:temperature[K],W:mixing ratio[kg/kg],P:pressure[Pa])