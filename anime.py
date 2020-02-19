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
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import os

######どのアニメを出力するか#####
# 0:微小角 1:カオス
CONFIG = 1
###############################

if(CONFIG == 0):
    FILE_NAME = "small_angle_anime"
if(CONFIG == 1):
    FILE_NAME = "chaos_anime"


df1 = pd.read_csv("data/" + FILE_NAME + ".csv")

t = df1["t"]
x1 = df1["x1"]
y1 = df1["y1"]
x2 = df1["x2"]
y2 = df1["y2"]


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
#ax.legend(fontsize=20)
# ax2.set_yscale("log")

##########################################################################
#animate 関数の中でこいつらだけを更新する
line, = ax.plot([], [], "o-", linewidth=2)
locus, = ax.plot([], [], "r-", linewidth=2)
time_text = ax.text(0.05, 0.9, "", transform=ax.transAxes, fontsize="20")
##########################################################################

xlocus, ylocus = [], []

def animate(cnt):
    x = [0, x1[cnt], x2[cnt]]
    y = [0, y1[cnt], y2[cnt]]
    line.set_data(x, y)

    xlocus.append(x2[cnt])
    ylocus.append(y2[cnt])
    locus.set_data(xlocus, ylocus)
  
    print(cnt) #確認
    time_text.set_text("time = %.2fs" % t[cnt])
    # return はblit=Trueの時のみ必要。Falseの時にあっても問題ない
    return line, locus, time_text


ani = FuncAnimation(fig, animate, frames=int(len(t)), interval=20, repeat=False, blit=True)

plt.show()
#ani.save("animes/anime.gif", writer="pillow", fps=20)
