import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def DrawArchimedean(a,b,c,phi,t0,t1,ax,v_num=1000):
    # 极角 t 的取值范围
    t = np.linspace(t0, t1, v_num)  # 多转几圈，看曲线效果

    # 极径 r 计算
    r = a + b * (t + phi)**(1/c)

    # 转换为笛卡尔坐标用于绘图
    x = r * np.cos(t)
    y = r * np.sin(t)

    # 颜色渐变，利用 hsv 色相值从0到1变化
    hues = np.linspace(0, 1, len(t))
    colors = [mcolors.hsv_to_rgb((h, 1, 1)) for h in hues]

    # 绘制曲线，逐点绘制颜色渐变

    for i in range(len(t)-1):
        ax.plot(x[i:i+2], y[i:i+2], color=colors[i], linewidth=2)
    

if __name__ == '__main__':
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    DrawArchimedean(1.0,2.0,-2.0,0.0,np.pi/2,np.pi*4,ax)
    plt.show()