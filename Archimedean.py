import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import random


def DrawArchimedean(a,b,c,phi,t0,t1,ax,rotate = 0.0,trans = np.array([0.0,0.0]),v_num=1000):
    # 极角 t 的取值范围
    t = np.linspace(t0, t1, v_num)  # 多转几圈，看曲线效果

    # 极径 r 计算
    r = a + b * (t + phi)**(1/c)

    # 转换为笛卡尔坐标用于绘图
    x = r * np.cos(t)
    y = r * np.sin(t)

    [x,y] = [x*np.cos(rotate)-y*np.sin(rotate),x*np.sin(rotate)+y*np.cos(rotate)]
    x = x+trans[0]
    y = y+trans[1]

    # 颜色渐变，利用 hsv 色相值从0到1变化
    hues = np.linspace(0, 1, len(t))
    colors = [mcolors.hsv_to_rgb((h, 1, 1)) for h in hues]

    # 绘制曲线，逐点绘制颜色渐变

    for i in range(len(t)-1):
        ax.plot(x[i:i+2], y[i:i+2], color=colors[i], linewidth=2)

def SignAngle(v1,v2):
    #计算从v1到v2的符号角
    product = np.dot(v1,v2)/np.linalg.norm(v1)/np.linalg.norm(v2)
    angle = np.arccos(product)
    if(v1[0]*v2[1]-v1[1]*v2[0]<0):
        angle = -angle
    return angle

def ComputeAngle(arc_length,a,b,c,phi):
    t0 = np.pi/2
    t1 = np.pi*4
    t = np.linspace(t0,t1,1000)
    alpha = []
    beta = []
    for st in t:
        et = st+arc_length
        sr = a + b * (st + phi)**(1/c)
        er = a + b * (et + phi)**(1/c)
        sx = sr * np.cos(st)
        sy = sr * np.sin(st)
        ex = er * np.cos(et)
        ey = er * np.sin(et)
        T = np.array([ex-sx,ey-sy])
        sr_d = b/c*(st+phi)**(1/c-1)
        er_d = b/c*(et+phi)**(1/c-1)
        s_tan = sr_d*np.array([np.cos(st),np.sin(st)])+sr*np.array([-np.sin(st),np.cos(st)])
        e_tan = er_d*np.array([np.cos(et),np.sin(et)])+er*np.array([-np.sin(et),np.cos(et)])
        alpha.append(SignAngle(s_tan,T)/np.pi*180.0)
        beta.append(SignAngle(T,e_tan)/np.pi*180.0)
    fig, ax = plt.subplots()
    ax.plot(t,alpha,color="red")
    ax.plot(t,beta,color="blue")

def GA_Inter_G0(ps,pe,c,center,ax):
    v1 = ps-center
    v2 = pe-center
    angle = SignAngle(v1,v2)
    ratio = (np.linalg.norm(v2)/np.linalg.norm(v1))**c
    t1 = angle/(ratio-1)
    t2 = angle*ratio/(ratio-1)
    if(t1<=0 or t2 <=0):
        return
    b = np.linalg.norm(v1)/(t1**(1/c))
    or_ps = (b*(t1**(1/c)))*np.array([np.cos(t1),np.sin(t1)])
    or_pe = (b*(t2**(1/c)))*np.array([np.cos(t2),np.sin(t2)])
    rotate_angle = SignAngle(or_pe-or_ps,pe-ps)
    x = or_ps[0]
    y = or_ps[1]
    new_or_ps= [x*np.cos(rotate_angle)-y*np.sin(rotate_angle),x*np.sin(rotate_angle)+y*np.cos(rotate_angle)]
    trans = ps-new_or_ps
    print(b,t1,t2)
    DrawArchimedean(0.0,b,c,0.0,t1,t2,ax,rotate_angle,trans)


if __name__ == '__main__':
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ps = np.array([-5.0,0.0])
    pe = np.array([12.0,3.0])
    ax.plot([ps[0],pe[0]],[ps[1],pe[1]],color="black")
    for i in range(5):
        random_x = random.uniform(0,4)
        random_y = random.uniform(-4,7)
        center = np.array([random_x,random_y])
        GA_Inter_G0(ps,pe,-2,center,ax)
    #DrawArchimedean(0.0,2.0,-2.0,0.0,np.pi/2,np.pi*7,ax)
    plt.show()