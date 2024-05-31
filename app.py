from flask import render_template,request,Flask
#from Steam import steam_calculator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as math
#import mpld3
import base64
import io

from scipy.interpolate import griddata

df_t= pd.read_csv("SST.csv")
df_p= pd.read_csv("SSP.csv")
df_c1= pd.read_csv("CLS.csv")
df_t.rename(columns={'T (C)':'Temperature'},inplace=True)
df_c=df_c1[df_c1["Phase"].isin(['Liquid','Vapor','SupercriticalFluid'])]
de=df_c1[df_c1["Phase"].isin(['SaturatedVapor','SaturatedLiquid'])]

df_c['Phase']=df_c['Phase'].astype(str)
df_c.dropna(inplace=True)
df_c=df_c.reset_index(drop=True)
df=df_c

dt = pd.read_csv("Entdry.csv")

def enth_dry(h,x):
    point = np.array(list(zip(dt['Enthalpy'],dt['dryness'])))
    value = np.array(list(zip(dt['Entropy_L'],dt['Entropy_V'],dt['Entropy'],dt['Enthalpy_L'],dt['Enthalpy_V'],)))
    q1 = griddata(point, value, (h,x), method='linear')
    ETRL=round(q1[0],3)
    ETRV=round(q1[1],3)
    ETR=round(q1[2],3)
    ETL=round(q1[3],3)
    ETV=round(q1[4],3)

    value = np.array(list(zip(df_t['Temperature'],df_t['P (MPa)'])))
    point = np.array(list(zip(df_t['Enthalpy Vapor (kJ/kg)'],df_t['Enthalpy Liquid (kJ/kg)'])))
    interpolated_value = griddata(point, value, (ETV,ETL), method='linear')
    k=interpolated_value
    t=round(k[0],3)
    p=round(k[1],3)
        
    plt.figure(figsize=(8, 6))
    
    
    plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=1)
    plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=1)
    
    # Set axis labels
    plt.xlabel("Entropy", color='black')
    plt.ylabel("Enthalpy", color='black')

    # Add a legend
    

    plt.scatter(ETR, h, c='black', marker='*', label='Specific Point', s=100, zorder =2)
    plt.grid(color='black', linestyle='--', linewidth=2)
    plt.legend()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
    image_data = base64.b64encode(buffer.getvalue()).decode()
    plt.close()
    return t,p,ETRL,ETRV,ETR,ETL,ETV



def entr_dry(s,x):
    point = np.array(list(zip(dt['Entropy'],dt['dryness'])))
    value = np.array(list(zip(dt['Entropy_L'],dt['Entropy_V'],dt['Enthalpy'],dt['Enthalpy_L'],dt['Enthalpy_V'],)))
    q1 = griddata(point, value, (s,x), method='linear')
    ETRL=round(q1[0],3)
    ETRV=round(q1[1],3)
    ET=round(q1[2],3)
    ETL=round(q1[3],3)
    ETV=round(q1[4],3)

    value = np.array(list(zip(df_t['Temperature'],df_t['P (MPa)'])))
    point = np.array(list(zip(df_t['Enthalpy Vapor (kJ/kg)'],df_t['Enthalpy Liquid (kJ/kg)'])))
    interpolated_value = griddata(point, value, (ETV,ETL), method='linear')
    k=interpolated_value
    t=round(k[0],3)
    p=round(k[1],3)
    plt.figure(figsize=(8, 6))
    
    
    plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=1)
    plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=1)
    
    # Set axis labels
    plt.xlabel("Entropy", color='black')
    plt.ylabel("Enthalpy", color='black')

    # Add a legend
    

    plt.scatter(s,ET, c='black', marker='*', label='Specific Point', s=100, zorder =2)
    plt.grid(color='black', linestyle='--', linewidth=2)
    plt.legend()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
    image_data = base64.b64encode(buffer.getvalue()).decode()
    plt.close()
    return t,p,ETRL,ETRV,ET,ETL,ETV



def steam_calculator(p,t):

    point = np.array(list(zip(de['Pressure(MPa)'])))
    value = np.array(list(zip(de['Temperature'],de['Specific_Enthalpy'])))
    q1 = griddata(point, value, (p), method='linear')
    Sat_Temp=round(float(q1[0]),2)


    if p<22.12 and t >= Sat_Temp :
        if t==Sat_Temp :
            Phase='Saturation Temperature'
        else :
            Phase='Vapor'

        Phase=Phase
        df=df_c1[df_c1["Phase"].isin(['Vapor','SaturatedVapor','SupercriticalFluid'])]
        points = np.array(list(zip(df['Pressure(MPa)'], df['Temperature'])))
        values = np.array(list(zip(df['Specific_Enthalpy'],df['Specific_Entropy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
        interpolated_value = griddata(points, values, (p,t), method='linear')
        k=interpolated_value

        Enthalpy=round(float(k[0]),2)
        Entropy=round(float(k[1]),2)
        Internal_energy=round(float(k[2]),2)
        Specific_volume=round(float(k[3]),5)
        Density=round(float(k[4]),2)

    elif p<22.12 and t <= Sat_Temp :
        if t==Sat_Temp :
            Phase='Saturation Temperature'
        else :
            Phase='Liquid'

        Phase=Phase
        df=df_c1[df_c1["Phase"].isin(['Liquid','SaturatedLiquid'])]
        points = np.array(list(zip(df['Pressure(MPa)'], df['Temperature'])))
        values = np.array(list(zip(df['Specific_Enthalpy'],df['Specific_Entropy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
        interpolated_value = griddata(points, values, (p,t), method='linear')
        k=interpolated_value

        Enthalpy=round(float(k[0]),2)
        Entropy=round(float(k[1]),2)
        Internal_energy=round(float(k[2]),2)
        Specific_volume=round(float(k[3]),5)
        Density=round(float(k[4]),2)

    elif p > 22.12 and t > 374 :

        Sat_Temp='No Saturation Temperature'

        Phase='Vapor'

        df=df_c1[df_c1["Phase"].isin(['Vapor','SupercriticalFluid'])]
        points = np.array(list(zip(df['Pressure(MPa)'], df['Temperature'])))
        values = np.array(list(zip(df['Specific_Enthalpy'],df['Specific_Entropy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
        interpolated_value = griddata(points, values, (p,t), method='linear')
        k=interpolated_value



        Enthalpy=round(float(k[0]),2)
        Entropy=round(float(k[1]),2)
        Internal_energy=round(float(k[2]),2)
        Specific_volume=round(float(k[3]),5)
        Density=round(float(k[4]),2)

    elif p > 22.12 and t < 374:
        Phase='Liquid'
        Sat_Temp='No Saturation Temperature'

        df=df_c1[df_c1["Phase"].isin(['Liquid'])]
        points = np.array(list(zip(df['Pressure(MPa)'], df['Temperature'])))
        values = np.array(list(zip(df['Specific_Enthalpy'],df['Specific_Entropy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
        interpolated_value = griddata(points, values, (p,t), method='linear')
        k=interpolated_value

        Enthalpy=round(float(k[0]),2)
        Entropy=round(float(k[1]),2)
        Internal_energy=round(float(k[2]),2)
        Specific_volume=round(float(k[3]),5)
        Density=round(float(k[4]),2)


    Ph=Phase
    plt.figure(figsize=(8, 6))
    
    plt.grid(color='black', linestyle='--', linewidth=2, zorder=1)
    plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=2)
    plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=2)
    
    # Set axis labels
    plt.xlabel("Entropy", color='black')
    plt.ylabel("Enthalpy", color='black')

    # Add a legend
    

    plt.scatter(Entropy, Enthalpy, c='black', marker='*', label='Specific Point', s=300, zorder =3)
    
    plt.legend()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
    image_data = base64.b64encode(buffer.getvalue()).decode()
    plt.close()



    return Enthalpy,Entropy,Ph,Sat_Temp,Specific_volume,Density,Internal_energy, image_data

def steam_entropy(p,s):

    if p<22.12 :
        point = np.array(list(zip(df_p['P (MPa)'])))
        value = np.array(list(zip(df_p['Entropy Vapor [kJ/(kg K)]'],df_p['Entropy Liquid [kJ/(kg K)]'])))
        interpolated_value = griddata(point, value, (p), method='linear')
        k=interpolated_value
        ETRV=k[0]
        ETRL=k[1]
        if s < ETRV and s > ETRL:
            q='L22.12 In'
            points = np.array(list(zip(df_p['P (MPa)'])))
            values = np.array(list(zip(df_p['T (C)'],df_p['Specific Volume Liquid (m^3/kg)'],df_p['Specific Volume Vapor (m^3/kg)'],df_p['Internal Energy Liquid (kJ/kg)'],df_p['Internal Energy Vapor (kJ/kg)'],df_p['Internal Energy of Vaporization (kJ/kg)'],df_p['Enthalpy Liquid (kJ/kg)'],df_p['Enthalpy Vapor (kJ/kg)'],df_p['Enthalpy of Vaporization (kJ/kg)'],df_p['Entropy Liquid [kJ/(kg K)]'],df_p['Entropy Vapor [kJ/(kg K)]'],df_p['Entropy of Vaporization [kJ/(kg K)]'])))
            interpolated_value = griddata(points, values, (p), method='linear')
            k=interpolated_value
            Sat_Temp=k[0]
            t=round(Sat_Temp,2)
            SVL=k[1]
            SVV=k[2]
            IEL=k[3]
            IEV=k[4]
            EL=k[6]
            EV=k[7]
            ETRL=k[9]
            ETRV=k[10]
            x = round((s-ETRL)/(ETRV-ETRL),3)
            Enthalpy=round(EL+x*(EV-EL),2)
            Specific_volume=round(SVL+x*(SVV-SVL),5)
            Internal_energy=round(IEL+x*(IEV-IEL),2)
            Density=round((1/Specific_volume),3)
            Ph = 'Saturated'
        elif s>ETRV or s<ETRL:
            q='L22.12 away'
            points = np.array(list(zip(df['Pressure(MPa)'], df['Specific_Entropy'])))
            values = np.array(list(zip(df['Temperature'],df['Specific_Enthalpy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
            interpolated_value = griddata(points, values, (p,s), method='linear')
            k=interpolated_value
            Phase=[]
            t=round(float(k[0]),2)
            Enthalpy=round(float(k[1]),2)
            Internal_energy=round(float(k[2]),2)
            Specific_volume=round(float(k[3]),5)
            Density=round(float(k[4]),2)
            point = np.array(list(zip(de['Pressure(MPa)'])))
            value = np.array(list(zip(de['Temperature'])))
            q1 = griddata(point, value, (p), method='linear')
            Sat_Temp=round(float(q1[0]),2)

            if p<22.12 and t > Sat_Temp :
                x=1
                Phase='Vapor'
            elif p<22.12 and t < Sat_Temp :
                x=0
                Phase='Liquid'
            elif p<22.12 and t==Sat_Temp :
                Phase='Saturation Temperature'
            Ph=Phase
    else:
        q='g22.12'
        points = np.array(list(zip(df['Pressure(MPa)'], df['Specific_Entropy'])))
        values = np.array(list(zip(df['Temperature'],df['Specific_Enthalpy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
        interpolated_value = griddata(points, values, (p,s), method='linear')
        k=interpolated_value
        Phase=[]
        t=round(float(k[0]),2)
        Enthalpy=round(float(k[1]),2)
        Internal_energy=round(float(k[2]),2)
        Specific_volume=round(float(k[3]),5)
        Density=round(float(k[4]),2)
        point = np.array(list(zip(de['Pressure(MPa)'])))
        value = np.array(list(zip(de['Temperature'])))
        q1 = griddata(point, value, (p), method='linear')
        Sat_Temp=round(float(q1[0]),2)

        if p > 22.12 and t > 374 :
            Phase='Super Critical Fluid'
            x=1
        elif p > 22.12 and t < 374:
            x=0
            Phase='Liquid'
            Sat_Temp='No Saturation Temperature'
        Ph=Phase
    plt.figure(figsize=(8, 6))


    plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=1)
    plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=1)
    
    # Set axis labels
    plt.xlabel("Entropy", color='black')
    plt.ylabel("Enthalpy", color='black')

    # Add a legend
    

    plt.scatter(s, Enthalpy, c='black', marker='*', label='Specific Point', s=100, zorder =2)
    plt.grid(color='black', linestyle='--', linewidth=2)
    plt.legend()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
    image_data = base64.b64encode(buffer.getvalue()).decode()
    plt.close()
    return t,Enthalpy,Ph,Sat_Temp,Specific_volume,Density,Internal_energy,x, image_data

def steam_enthalpy(p,h):

    if p<22.12 :
        point = np.array(list(zip(df_p['P (MPa)'])))
        value = np.array(list(zip(df_p['Enthalpy Vapor (kJ/kg)'],df_p['Enthalpy Liquid (kJ/kg)'])))
        interpolated_value = griddata(point, value, (p), method='linear')
        k=interpolated_value
        EV=k[0]
        EL=k[1]
        if h < EV and h > EL:
            q='L22.12 In'
            points = np.array(list(zip(df_p['P (MPa)'])))
            values = np.array(list(zip(df_p['T (C)'],df_p['Specific Volume Liquid (m^3/kg)'],df_p['Specific Volume Vapor (m^3/kg)'],df_p['Internal Energy Liquid (kJ/kg)'],df_p['Internal Energy Vapor (kJ/kg)'],df_p['Internal Energy of Vaporization (kJ/kg)'],df_p['Enthalpy Liquid (kJ/kg)'],df_p['Enthalpy Vapor (kJ/kg)'],df_p['Enthalpy of Vaporization (kJ/kg)'],df_p['Entropy Liquid [kJ/(kg K)]'],df_p['Entropy Vapor [kJ/(kg K)]'],df_p['Entropy of Vaporization [kJ/(kg K)]'])))
            interpolated_value = griddata(points, values, (p), method='linear')
            k=interpolated_value
            Sat_Temp=k[0]
            t=round(Sat_Temp,2)
            SVL=k[1]
            SVV=k[2]
            IEL=k[3]
            IEV=k[4]
            EL=k[6]
            EV=k[7]
            ETRL=k[9]
            ETRV=k[10]
            x = round((h-EL)/(EV-EL),3)
            Entropy=round(ETRL+x*(ETRV-ETRL),2)
            Specific_volume=round(SVL+x*(SVV-SVL),5)
            Internal_energy=round(IEL+x*(IEV-IEL),2)
            Density=round((1/Specific_volume),3)
            Ph = 'Saturated'
        elif h>EV or h<EL:
            q='L22.12 away'
            points = np.array(list(zip(df['Pressure(MPa)'], df['Specific_Enthalpy'])))
            values = np.array(list(zip(df['Temperature'],df['Specific_Entropy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
            interpolated_value = griddata(points, values, (p,h), method='linear')
            k=interpolated_value
            Phase=[]
            t=round(float(k[0]),2)
            Entropy=round(float(k[1]),2)
            Internal_energy=round(float(k[2]),2)

            Density=round(float(k[4]),2)
            Specific_volume=round((1/Density),5)

            point = np.array(list(zip(de['Pressure(MPa)'])))
            value = np.array(list(zip(de['Temperature'])))
            q1 = griddata(point, value, (p), method='linear')
            Sat_Temp=round(float(q1[0]),2)

            if p<22.12 and t > Sat_Temp :
                x=1
                Phase='Vapor'
            elif p<22.12 and t < Sat_Temp :
                x=0
                Phase='Liquid'
            elif p<22.12 and t==Sat_Temp :
                Phase='Saturation Temperature'
            Ph=Phase
    else:
        q='g22.12'
        points = np.array(list(zip(df['Pressure(MPa)'], df['Specific_Enthalpy'])))
        values = np.array(list(zip(df['Temperature'],df['Specific_Entropy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
        interpolated_value = griddata(points, values, (p,h), method='linear')
        k=interpolated_value
        Phase=[]
        t=round(float(k[0]),2)
        Entropy=round(float(k[1]),2)
        Internal_energy=round(float(k[2]),2)

        Density=round(float(k[4]),2)
        Specific_volume=round((1/Density),5)
        point = np.array(list(zip(de['Pressure(MPa)'])))
        value = np.array(list(zip(de['Temperature'])))
        q1 = griddata(point, value, (p), method='linear')
        Sat_Temp=round(float(q1[0]),2)

        if p > 22.12 and t > 374 :
            x=1
            Phase='Super Critical Fluid'
        elif p > 22.12 and t < 374:
            x=0
            Phase='Liquid'
            Sat_Temp='No Saturation Temperature'
        Ph=Phase
        plt.figure(figsize=(8, 6))
        
        
        plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=1)
        plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=1)
        
        # Set axis labels
        plt.xlabel("Entropy", color='black')
        plt.ylabel("Enthalpy", color='black')

        # Add a legend
        

        plt.scatter(Entropy, h, c='black', marker='*', label='Specific Point', s=100, zorder =2)
        plt.grid(color='black', linestyle='--', linewidth=2)
        plt.legend()

        buffer = io.BytesIO()
        plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
        image_data = base64.b64encode(buffer.getvalue()).decode()
        plt.close()
    return t,Entropy,Ph,Sat_Temp,Specific_volume,Density,Internal_energy,x, image_data


def steam_dry(p,x):
    points = np.array(list(zip(df_p['P (MPa)'])))
    values = np.array(list(zip(df_p['T (C)'],df_p['Specific Volume Liquid (m^3/kg)'],df_p['Specific Volume Vapor (m^3/kg)'],df_p['Internal Energy Liquid (kJ/kg)'],df_p['Internal Energy Vapor (kJ/kg)'],df_p['Internal Energy of Vaporization (kJ/kg)'],df_p['Enthalpy Liquid (kJ/kg)'],df_p['Enthalpy Vapor (kJ/kg)'],df_p['Enthalpy of Vaporization (kJ/kg)'],df_p['Entropy Liquid [kJ/(kg K)]'],df_p['Entropy Vapor [kJ/(kg K)]'],df_p['Entropy of Vaporization [kJ/(kg K)]'])))
    interpolated_value = griddata(points, values, (p), method='linear')
    k=interpolated_value
    t=round(k[0],2)
    SVL=k[1]
    SVV=k[2]
    SV=round(SVL+x*(k[2]-k[1]),5)
    IEL=k[3]
    IEV=k[4]
    IE=round(IEL+x*(k[4]-k[3]),2)
    EL=k[6]
    EV=k[7]
    EVW=k[8]
    E=round(EL+x*(EV-EL),2)
    ETRL=k[9]
    ETRV=k[10]
    ETR=round(ETRL+x*(ETRV-ETRL),2)
    ETRVW=k[11]
    plt.figure(figsize=(8, 6))
    
    
    plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=1)
    plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=1)
    
    # Set axis labels
    plt.xlabel("Entropy", color='black')
    plt.ylabel("Enthalpy", color='black')

    # Add a legend
    

    plt.scatter(ETR, E, c='black', marker='*', label='Specific Point', s=100, zorder =2)
    plt.grid(color='black', linestyle='--', linewidth=2)
    plt.legend()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
    image_data = base64.b64encode(buffer.getvalue()).decode()
    plt.close()
    return t,E,SV,ETR,IE, image_data

def temp_enthalpy(t,h):

    if t<374 :
        point = np.array(list(zip(df_t['Temperature'])))
        value = np.array(list(zip(df_t['Enthalpy Vapor (kJ/kg)'],df_t['Enthalpy Liquid (kJ/kg)'])))
        interpolated_value = griddata(point, value, (t), method='linear')
        k=interpolated_value
        EV=k[0]
        EL=k[1]
        if h < EV and h > EL:
            q='L22.12 In'
            points = np.array(list(zip(df_t['Temperature'])))
            values = np.array(list(zip(df_t['P (MPa)'],df_t['Specific Volume Liquid (m^3/kg)'],df_t['Specific Volume Vapor (m^3/kg)'],df_t['Internal Energy Liquid (kJ/kg)'],df_t['Internal Energy Vapor (kJ/kg)'],df_t['Internal Energy of Vaporization (kJ/kg)'],df_t['Enthalpy Liquid (kJ/kg)'],df_t['Enthalpy Vapor (kJ/kg)'],df_t['Enthalpy of Vaporization'],df_t['Entropy Liquid [kJ/(kg K)]'],df_t['Entropy Vapor [kJ/(kg K)]'],df_t['Entropy of Vaporization [kJ/(kg K)]'])))
            interpolated_value = griddata(points, values, (t), method='linear')
            k=interpolated_value
            Sat_p=round(k[0],2)
            p=round(Sat_p,2)
            SVL=k[1]
            SVV=k[2]
            IEL=k[3]
            IEV=k[4]
            EL=k[6]
            EV=k[7]
            ETRL=k[9]
            ETRV=k[10]
            x = round((h-EL)/(EV-EL),3)
            Entropy=round(ETRL+x*(ETRV-ETRL),2)
            Specific_volume=round(SVL+x*(SVV-SVL),5)
            Internal_energy=round(IEL+x*(IEV-IEL),2)
            Density=round((1/Specific_volume),3)
            Ph = 'Saturated'
        elif h>EV or h<EL:
            q='L22.12 away'
            points = np.array(list(zip(df['Temperature'], df['Specific_Enthalpy'])))
            values = np.array(list(zip(df['Pressure(MPa)'],df['Specific_Entropy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
            interpolated_value = griddata(points, values, (t,h), method='linear')
            k=interpolated_value
            Phase=[]
            p=round(float(k[0]),2)
            Entropy=round(float(k[1]),2)
            Internal_energy=round(float(k[2]),2)

            Density=round(float(k[4]),2)
            Specific_volume=round((1/Density),5)

            point = np.array(list(zip(df_t['Temperature'])))
            value = np.array(list(zip(df_t['P (MPa)'])))
            q1 = griddata(point, value, (t), method='linear')
            Sat_p=round(float(q1[0]),2)

            if h > EV:
                x=1
                Phase='Vapor'
            elif h<EL :
                x=0
                Phase='Liquid'
            elif h==EV or h==EL:
                Phase='Saturation Temperature'
            Ph=Phase
    else:
        q='g22.12'
        points = np.array(list(zip(df['Temperature'], df['Specific_Enthalpy'])))
        values = np.array(list(zip(df['Pressure(MPa)'],df['Specific_Entropy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
        interpolated_value = griddata(points, values, (t,h), method='linear')
        k=interpolated_value
        Phase=[]
        p=round(float(k[0]),2)
        Entropy=round(float(k[1]),2)
        Internal_energy=round(float(k[2]),2)

        Density=round(float(k[4]),2)
        Specific_volume=round((1/Density),5)
        point = np.array(list(zip(df_t['Temperature'])))
        value = np.array(list(zip(df_t['P (MPa)'])))
        q1 = griddata(point, value, (t), method='linear')
        Sat_p=round(float(q1[0]),2)

        if p > 22.12 and t > 374 :
            x=1
            Phase='Super Critical Fluid'
        elif p > 22.12 and t < 374:
            x=0
            Phase='Liquid'
            Sat_Temp='No Saturation Temperature'
        Ph=Phase
        plt.figure(figsize=(8, 6))
    
    
        plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=1)
        plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=1)
        
        # Set axis labels
        plt.xlabel("Entropy", color='black')
        plt.ylabel("Enthalpy", color='black')

        # Add a legend
        

        plt.scatter(Entropy, h, c='black', marker='*', label='Specific Point', s=100, zorder =2)
        plt.grid(color='black', linestyle='--', linewidth=2)
        plt.legend()

        buffer = io.BytesIO()
        plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
        image_data = base64.b64encode(buffer.getvalue()).decode()
        plt.close()
    return p,Entropy,Ph,Sat_p,Specific_volume,Density,Internal_energy,x, image_data

def temp_entropy(t,s):

    if t<374 :
        point = np.array(list(zip(df_t['Temperature'])))
        value = np.array(list(zip(df_t['Entropy Vapor [kJ/(kg K)]'],df_t['Entropy Liquid [kJ/(kg K)]'])))
        interpolated_value = griddata(point, value, (t), method='linear')
        k=interpolated_value
        ETRV=k[0]
        ETRL=k[1]
        if s < ETRV and s > ETRL:
            q='L22.12 In'
            points = np.array(list(zip(df_t['Temperature'])))
            values = np.array(list(zip(df_t['P (MPa)'],df_t['Specific Volume Liquid (m^3/kg)'],df_t['Specific Volume Vapor (m^3/kg)'],df_t['Internal Energy Liquid (kJ/kg)'],df_t['Internal Energy Vapor (kJ/kg)'],df_t['Internal Energy of Vaporization (kJ/kg)'],df_t['Enthalpy Liquid (kJ/kg)'],df_t['Enthalpy Vapor (kJ/kg)'],df_t['Enthalpy of Vaporization'],df_t['Entropy Liquid [kJ/(kg K)]'],df_t['Entropy Vapor [kJ/(kg K)]'],df_t['Entropy of Vaporization [kJ/(kg K)]'])))
            interpolated_value = griddata(points, values, (t), method='linear')
            k=interpolated_value
            Sat_p=round(k[0],2)
            p=round(Sat_p,2)
            SVL=k[1]
            SVV=k[2]
            IEL=k[3]
            IEV=k[4]
            EL=k[6]
            EV=k[7]
            ETRL=k[9]
            ETRV=k[10]
            x = round((s-ETRL)/(ETRV-ETRL),3)
            Enthalpy=round(EL+x*(EV-EL),2)
            Specific_volume=round(SVL+x*(SVV-SVL),5)
            Internal_energy=round(IEL+x*(IEV-IEL),2)
            Density=round((1/Specific_volume),3)
            Ph = 'Saturated'
        elif s>ETRV or s<ETRL:
            q='L22.12 away'
            points = np.array(list(zip(df['Temperature'], df['Specific_Entropy'])))
            values = np.array(list(zip(df['Pressure(MPa)'],df['Specific_Enthalpy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
            interpolated_value = griddata(points, values, (t,s), method='linear')
            k=interpolated_value
            Phase=[]
            p=round(float(k[0]),2)
            Enthalpy=round(float(k[1]),2)
            Internal_energy=round(float(k[2]),2)

            Density=round(float(k[4]),2)
            Specific_volume=round((1/Density),5)

            point = np.array(list(zip(df_t['Temperature'])))
            value = np.array(list(zip(df_t['P (MPa)'])))
            q1 = griddata(point, value, (t), method='linear')
            Sat_p=round(float(q1[0]),2)

            if s > ETRV:
                x=1
                Phase='Vapor'
            elif s < ETRL :
                x=0
                Phase='Liquid'
            elif s==ETRV or s==ETRL:
                Phase='Saturation Temperature'
            Ph=Phase
    else:
        q='g22.12'
        points = np.array(list(zip(df['Temperature'], df['Specific_Entropy'])))
        values = np.array(list(zip(df['Pressure(MPa)'],df['Specific_Enthalpy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
        interpolated_value = griddata(points, values, (t,s), method='linear')
        k=interpolated_value
        Phase=[]
        p=round(float(k[0]),2)
        Enthalpy=round(float(k[1]),2)
        Internal_energy=round(float(k[2]),2)

        Density=round(float(k[4]),2)
        Specific_volume=round((1/Density),5)
        point = np.array(list(zip(df_t['Temperature'])))
        value = np.array(list(zip(df_t['P (MPa)'])))
        q1 = griddata(point, value, (t), method='linear')
        Sat_p=round(float(q1[0]),2)

        if p > 22.12 and t > 374 :
            x=1
            Phase='Super Critical Fluid'
        elif p > 22.12 and t < 374:
            x=0
            Phase='Liquid'
            Sat_Temp='No Saturation Temperature'
        Ph=Phase
        plt.figure(figsize=(8, 6))
    
    
        plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=1)
        plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=1)
        
        # Set axis labels
        plt.xlabel("Entropy", color='black')
        plt.ylabel("Enthalpy", color='black')

        # Add a legend
        

        plt.scatter(s, Enthalpy, c='black', marker='*', label='Specific Point', s=100, zorder =2)
        plt.grid(color='black', linestyle='--', linewidth=2)
        plt.legend()

        buffer = io.BytesIO()
        plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
        image_data = base64.b64encode(buffer.getvalue()).decode()
        plt.close()
    return p,Enthalpy,Ph,Sat_p,Specific_volume,Density,Internal_energy,x, image_data

def temp_dry(t,x):
    q='L22.12 In'
    points = np.array(list(zip(df_t['Temperature'])))
    values = np.array(list(zip(df_t['P (MPa)'],df_t['Specific Volume Liquid (m^3/kg)'],df_t['Specific Volume Vapor (m^3/kg)'],df_t['Internal Energy Liquid (kJ/kg)'],df_t['Internal Energy Vapor (kJ/kg)'],df_t['Internal Energy of Vaporization (kJ/kg)'],df_t['Enthalpy Liquid (kJ/kg)'],df_t['Enthalpy Vapor (kJ/kg)'],df_t['Enthalpy of Vaporization'],df_t['Entropy Liquid [kJ/(kg K)]'],df_t['Entropy Vapor [kJ/(kg K)]'],df_t['Entropy of Vaporization [kJ/(kg K)]'])))
    interpolated_value = griddata(points, values, (t), method='linear')
    k=interpolated_value
    Sat_p=k[0]
    p=round(Sat_p,2)
    SVL=k[1]
    SVV=k[2]
    IEL=k[3]
    IEV=k[4]
    EL=k[6]
    EV=k[7]
    ETRL=k[9]
    ETRV=k[10]
    Entropy = round(ETRL+x*(ETRV-ETRL),2)
    Enthalpy=round(EL+x*(EV-EL),2)
    Specific_volume=round(SVL+x*(SVV-SVL),5)
    Internal_energy=round(IEL+x*(IEV-IEL),2)
    Density=round((1/Specific_volume),3)
    Ph = 'Saturated'
    
    plt.figure(figsize=(8, 6))
    
    
    plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=1)
    plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=1)
    
    # Set axis labels
    plt.xlabel("Entropy", color='black')
    plt.ylabel("Enthalpy", color='black')

    # Add a legend
    

    plt.scatter(Entropy, Enthalpy, c='black', marker='*', label='Specific Point', s=100, zorder =2)
    plt.grid(color='black', linestyle='--', linewidth=2)
    plt.legend()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
    image_data = base64.b64encode(buffer.getvalue()).decode()
    plt.close()
    return Sat_p,Enthalpy,Specific_volume,Entropy,Internal_energy, image_data


def enth_entropy(h,s):
    point = np.array(list(zip(dt['Enthalpy'],dt['Entropy'])))
    value = np.array(list(zip(dt['dryness'],dt['Enthalpy_L'],dt['Enthalpy_V'],dt['Entropy_L'],dt['Entropy_V'])))
    q1 = griddata(point, value, (h,s), method='linear')
    x1=q1[0]
    ETL=q1[1]
    ETV=q1[2]
    ETRL=q1[3]
    ETRV=q1[4]


    if math.isnan(x1):
        #print("The value does not exist or is NaN.")
        points = np.array(list(zip(df['Specific_Enthalpy'], df['Specific_Entropy'])))
        values = np.array(list(zip(df['Pressure(MPa)'],df['Temperature'])))
        interpolated_value = griddata(points, values, (h,s), method='linear')
        k=interpolated_value
        Phase=[]
        p=round(float(k[0]),2)
        t=round(float(k[1]),2)
        if p<22.12 :
            point = np.array(list(zip(df_p['P (MPa)'])))
            value = np.array(list(zip(df_p['Entropy Vapor [kJ/(kg K)]'],df_p['Entropy Liquid [kJ/(kg K)]'])))
            interpolated_value = griddata(point, value, (p), method='linear')
            k=interpolated_value
            ETRV=k[0]
            ETRL=k[1]
            points = np.array(list(zip(df['Pressure(MPa)'], df['Specific_Entropy'])))
            values = np.array(list(zip(df['Temperature'],df['Specific_Enthalpy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
            interpolated_value = griddata(points, values, (p,s), method='linear')
            k=interpolated_value
            Phase=[]
            t=round(float(k[0]),2)
            Enthalpy=round(float(k[1]),2)
            Internal_energy=round(float(k[2]),2)
            Specific_volume=round(float(k[3]),5)
            Density=round(float(k[4]),2)
            point = np.array(list(zip(de['Pressure(MPa)'])))
            value = np.array(list(zip(de['Temperature'])))
            q1 = griddata(point, value, (p), method='linear')
            Sat_Temp=round(float(q1[0]),2)

            if p<22.12 and t > Sat_Temp :
                x=1
                Phase='Vapor'
            elif p<22.12 and t < Sat_Temp :
                x=0
                Phase='Liquid'
            elif p<22.12 and t==Sat_Temp :
                Phase='Saturation Temperature'
            Ph=Phase
        else:
            q='g22.12'
            points = np.array(list(zip(df['Pressure(MPa)'], df['Specific_Entropy'])))
            values = np.array(list(zip(df['Temperature'],df['Specific_Enthalpy'],df['Specific_Internal_Energy'],df['Specific_Volume'],df['Density'])))
            interpolated_value = griddata(points, values, (p,s), method='linear')
            k=interpolated_value
            Phase=[]
            t=round(float(k[0]),2)
            Enthalpy=round(float(k[1]),2)
            Internal_energy=round(float(k[2]),2)
            Specific_volume=round(float(k[3]),5)
            Density=round(float(k[4]),2)
            point = np.array(list(zip(de['Pressure(MPa)'])))
            value = np.array(list(zip(de['Temperature'])))
            q1 = griddata(point, value, (p), method='linear')
            Sat_Temp=round(float(q1[0]),2)

            if p > 22.12 and t > 374 :
                Phase='Super Critical Fluid'
                Sat_Temp='No Saturation Temperature'
                x=1
            elif p > 22.12 and t < 374:
                x=0
                Phase='Liquid'
                Sat_Temp='No Saturation Temperature'
            Ph=Phase


    else:

        point = np.array(list(zip(dt['Enthalpy'],dt['Entropy'])))
        value = np.array(list(zip(dt['dryness'],dt['Enthalpy_L'],dt['Enthalpy_V'],dt['Entropy_L'],dt['Entropy_V'])))
        q1 = griddata(point, value, (h,s), method='linear')
        x=round(q1[0],3)
        ETL=q1[1]
        ETV=q1[2]
        ETRL=q1[3]
        ETRV=q1[4]

        value = np.array(list(zip(df_t['Temperature'],df_t['P (MPa)'])))
        point = np.array(list(zip(df_t['Enthalpy Vapor (kJ/kg)'],df_t['Enthalpy Liquid (kJ/kg)'])))
        interpolated_value = griddata(point, value, (ETV,ETL), method='linear')
        k=interpolated_value
        t=round(k[0],3)
        p=round(k[1],3)

        Ph= 'Saturated'
        points = np.array(list(zip(df_p['P (MPa)'],df_p['T (C)'])))
        values = np.array(list(zip(df_p['Specific Volume Liquid (m^3/kg)'],df_p['Specific Volume Vapor (m^3/kg)'],df_p['Internal Energy Liquid (kJ/kg)'],df_p['Internal Energy Vapor (kJ/kg)'],df_p['Internal Energy of Vaporization (kJ/kg)'],df_p['Enthalpy Liquid (kJ/kg)'],df_p['Enthalpy Vapor (kJ/kg)'],df_p['Enthalpy of Vaporization (kJ/kg)'],df_p['Entropy Liquid [kJ/(kg K)]'],df_p['Entropy Vapor [kJ/(kg K)]'],df_p['Entropy of Vaporization [kJ/(kg K)]'])))
        interpolated_value = griddata(points, values, (p,t), method='linear')
        k=interpolated_value
        Sat_Temp=k[0]
        #t=round(Sat_Temp,2)
        SVL=k[0]
        SVV=k[1]
        IEL=k[2]
        IEV=k[3]
        EL=k[5]
        EV=k[6]

        Enthalpy=round(EL+x*(EV-EL),2)
        Specific_volume=round(SVL+x*(SVV-SVL),5)
        Internal_energy=round(IEL+x*(IEV-IEL),2)
        Density=round((1/Specific_volume),3)
    
    plt.figure(figsize=(8, 6))
    
    
    plt.plot(df_t['Entropy Vapor [kJ/(kg K)]'], df_t['Enthalpy Vapor (kJ/kg)'], 'ro-', label='Vapor line', markersize=10, zorder=1)
    plt.plot(df_t['Entropy Liquid [kJ/(kg K)]'], df_t['Enthalpy Liquid (kJ/kg)'], 'bo-', label='Liquid line', markersize=10, zorder=1)
    
    # Set axis labels
    plt.xlabel("Entropy", color='black')
    plt.ylabel("Enthalpy", color='black')

    # Add a legend
    

    plt.scatter(s, h, c='black', marker='*', label='Specific Point', s=100, zorder =2)
    plt.grid(color='black', linestyle='--', linewidth=2)
    plt.legend()

    buffer = io.BytesIO()
    plt.savefig(buffer, format="png", dpi=300)  # Adjust DPI as needed
    image_data = base64.b64encode(buffer.getvalue()).decode()
    plt.close()

    return p,t,Ph,Sat_Temp,Specific_volume,Density,Internal_energy,x, image_data




app=Flask(__name__)









@app.route('/')
@app.route('/home')
def home():
    return render_template("PT.html")

@app.route('/result', methods=['GET','POST'])

def index():

    result, result_2, result_3, result_4, result_5, result_6, result_7, result_8, result_9, result_10 = [None] * 10

    if request.method=='POST':
        function_name=request.form['condition']
		

        if function_name=='PT':
            pressure_unit = request.form['units_pressure']
            temperature_unit =request.form['units_temperature']
            p=float(request.form['pressure'])
            t=float(request.form['temperature'])
            if pressure_unit == "BAR":
                p= p*0.1

            if temperature_unit == "F":
                t = (t - 32) * 5/9

            result=steam_calculator(p,t)


        elif function_name=='PE':
            pressure_unit1 = request.form['units_pressure1']
            entropy_unit =request.form['units_entropy']
            p=float(request.form['pressure1'])
            s=float(request.form['entropy'])
            if pressure_unit1 == "BAR":
                p= p*0.1

            if entropy_unit == "J/kg-K":
                s= s/1000
            result_2=steam_entropy(p,s)


        elif function_name=='PH':
            pressure_unit2 = request.form['units_pressure2']
            enthalpy_unit =request.form['units_enthalpy']
            p=float(request.form['pressure2'])
            h=float(request.form['enthalpy'])
            if pressure_unit2 == "BAR":
                p= p*0.1
            if enthalpy_unit == "J/kg":
                h=h/1000
            result_3=steam_enthalpy(p,h)

        elif function_name=='PX':
            pressure_unit3 = request.form['units_pressure3']
            dryness_unit =request.form['units_dryness']
            p=float(request.form['pressure3'])
            x=float(request.form['dryness'])
            if pressure_unit3 == "BAR":
                p= p*0.1
            if dryness_unit == "%":
                x = x/100
            result_4=steam_dry(p,x)
        elif function_name=='TE':
            temperature_unit1 = request.form['units_temperature1']
            entropy_unit1 =request.form['units_entropy1']
            t=float(request.form['temperature1'])
            s=float(request.form['entropy1'])
            if temperature_unit1 == "F":
                t = (t - 32) * 5/9
            if entropy_unit1 == "J/kg-K":
                s = s/1000
            result_5= temp_entropy(t,s)
        elif function_name=='TH':
            temperature_unit2 = request.form['units_temperature2']
            enthalpy_unit1 =request.form['units_enthalpy1']
            t=float(request.form['temperature2'])
            h=float(request.form['enthalpy1'])
            if temperature_unit2 == "F":
                t = (t - 32) * 5/9
            if enthalpy_unit1 == "J/kg":
                h = h/1000
            result_6= temp_enthalpy(t,h)

        elif function_name=='TX':
            temperature_unit3 = request.form['units_temperature3']
            dryness_unit1 =request.form['units_dryness1']
            t=float(request.form['temperature3'])
            x=float(request.form['dryness1'])
            if temperature_unit3 == "F":
                t = (t - 32) * 5/9
            if dryness_unit1 == "%":
                x = x/100
            result_7= temp_dry(t,x)

        elif function_name=='HE':
            enthalpy_unit2 =request.form['units_enthalpy2']
            entropy_unit2 =request.form['units_entropy2']
            h= float(request.form['enthalpy2'])
            s= float(request.form['entropy2'])
            if enthalpy_unit2 == "J/kg":
                h = h/1000
            if entropy_unit2 == "J/kg-K":
                s = s/1000
            result_8 = enth_entropy(h,s)

        elif function_name=='HX':
            enthalpy_unit3 =request.form['units_enthalpy3']
            dryness_unit2 =request.form['units_dryness2']
            h= float(request.form['enthalpy3'])
            x= float(request.form['dryness2'])
            if enthalpy_unit3 == "J/kg":
                h = h/1000
            if dryness_unit2 == "%":
                x = x/100
            result_9 = enth_dry(h,x)

        elif function_name=='EX':
            entropy_unit3 =request.form['units_entropy3']
            dryness_unit3 =request.form['units_dryness3']
            s= float(request.form['entropy3'])
            x= float(request.form['dryness3'])
            if entropy_unit3 == "J/kg-K":
                s = s/1000
            if dryness_unit3 == "%":
                x = x/100
            result_10 = entr_dry(s,x)




    return render_template("PT.html", result=result, result_2=result_2, result_3=result_3, result_4=result_4,result_5=result_5,result_6=result_6,result_7=result_7,result_8=result_8,result_9=result_9,result_10=result_10)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=3000)
