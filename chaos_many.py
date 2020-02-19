#############使い方###################
'''
相対パス "data/chaos_anime.csv" の中に
t,x1,y1,x2,y2
(以下データ)

の形で位置情報を入れておけばこれを動かせばアニメが出来る

[注意]
GIFだと写真の枚数600枚くらいが限界
⇒データの出力は適度な頻度に抑える必要がある
例
写真の枚数を500枚で20[s]の軌跡を追うために
20/500[s]ごとのデータをcsvファイルに出力する
'''
######################################

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

import os

############################CONFIG############################
GIF = 1
MP4 = 0
PLT = 0
##############################################################


######初期条件の個数#####
nfp = open("data_many/init_num.txt")
NUM = int(nfp.readline())
nfp.close()
###############################



FILE_NAME = "many_chaos"

df = pd.read_csv("data_many/" + FILE_NAME + "0.csv")
t=df["t"]
x1=[]
y1=[]
x2=[]
y2=[]

for i in range(NUM):
    df = pd.read_csv("data_many/" + FILE_NAME + str(i) + ".csv")
    x1.append(df["x1"])
    y1.append(df["y1"])
    x2.append(df["x2"])
    y2.append(df["y2"])

fig = plt.figure(figsize=(8, 8))
# fig_xy.subplots_adjust(left=0.80)
ax = fig.add_subplot(111)
#ax2.hlines([0], tmin, tmax)
ax.set_ylabel("y", fontsize=18)
ax.set_xlabel("x", fontsize=18)
ax.grid(linestyle="dotted")
ax.xaxis.set_tick_params(direction='in')
ax.yaxis.set_tick_params(direction='in')
ax.set_xlim(-2.1, 2.1)
ax.set_ylim(-2.1, 2.1)
ax.tick_params(labelsize=14)  # 軸目盛の数字のサイズ
ax.set_title("Double Pendulums", fontsize=20)
#ax.legend(fontsize=20)
# ax2.set_yscale("log")


#line, = ax.plot([], [], "o-", linewidth=2)

##########################################################################
#animate 関数の中でこいつらだけを更新する

lines=[]
locus=[]
for i in range(NUM):
    line, = ax.plot([], [], "o-", linewidth=2, color=cm.hsv(i/NUM))
    loc, = ax.plot([], [], linewidth=2)
    lines.append(line)
    locus.append(loc)

time_text = ax.text(0.05, 0.9, "", transform=ax.transAxes, fontsize="20")
##########################################################################



xlocus, ylocus = [], []
for _ in range(NUM):
    xlocus.append([])
    ylocus.append([])

def animate(cnt):
    for i in range(NUM):
        #xlocus[i].append(x2[i][cnt])
        #ylocus[i].append(y2[i][cnt])
        #locus[i].set_data(xlocus[i], ylocus[i])

        x = [0, x1[i][cnt], x2[i][cnt]]
        y = [0, y1[i][cnt], y2[i][cnt]]
        lines[i].set_data(x, y)
    

    print(cnt) #確認
    time_text.set_text("time = %.2fs" % t[cnt])


ani = FuncAnimation(fig, animate, frames=int(len(t)), interval=20, repeat=True, blit=False)

if(GIF == 1):
    ani.save("animes/anime.gif", writer="pillow", fps=20)
if(MP4 == 1):
    ani.save("animes/anime.mp4", writer="ffmpeg", fps=20)
if(PLT == 1):
    plt.show()
