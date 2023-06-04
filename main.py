# -*- coding: cp1251 -*-
import tkinter as tk
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def create():
    global swarm
    global velocity
    global pbest_pos
    global pbest_value
    global gbest_pos
    global gbest_value
    global k
    global bounds
    global dim
    global mystat
    mystat = []
    ax.clear()
    canvas.draw()
    dim = int(mutant_entry.get())
    k_entry.delete(0, tk.END)
    k_entry.insert(0, 0)
    k = 0
    bounds = np.array([int(minn_entry.get()) * dim, int(maxx_entry.get()) * dim])
    swarm = np.random.uniform(bounds[0], bounds[1], (int(chromoval_entry.get()), dim))
    velocity = np.zeros((int(chromoval_entry.get()), dim))
    pbest_pos = swarm.copy()
    pbest_value = np.array([np.inf] * int(chromoval_entry.get()))
    gbest_pos = np.zeros(dim)
    gbest_value = np.inf

def draw_graph(points):
    ax.clear()
    x = [point[0] for point in points]
    y = [point[1] for point in points]
    ax.scatter(x, y)
    canvas.draw()

def function(x):
    return 8*x[0]**2+4*x[0]*x[1]+5*x[1]**2

def main():
    global swarm
    global velocity
    global pbest_pos
    global pbest_value
    global gbest_pos
    global gbest_value
    global k
    global bounds
    global mystat
    gener_val = int(generations_var.get())
    chromo_val = int(chromoval_entry.get())
    k += gener_val
    k_entry.delete(0, tk.END)
    k_entry.insert(0, k)
    w = float(inercia_entry.get())
    c1 = c2 = float(speed_entry.get())
    for i in range(gener_val):
        f = np.array([function(x) for x in swarm])

        mask = f < pbest_value
        pbest_value[mask] = f[mask]
        pbest_pos[mask] = swarm[mask]

        mask = pbest_value < gbest_value
        if np.any(mask):
            gbest_value = np.min(pbest_value)
            gbest_pos = pbest_pos[np.argmin(pbest_value)]
        if gbest_value < 80 :
            mystat.append(gbest_value)

        r1 = np.random.rand(chromo_val, dim)
        r2 = np.random.rand(chromo_val, dim)
        velocity = w * velocity + c1 * r1 * (pbest_pos - swarm) + c2 * r2 * (gbest_pos - swarm)
        swarm = swarm + velocity

        swarm = np.clip(swarm, bounds[0], bounds[1])

        swarm = np.where(swarm < bounds[0], bounds[0], swarm)
        swarm = np.where(swarm > bounds[1], bounds[1], swarm)
    draw_graph(swarm)
    output_text2.delete(1.0, tk.END)
    output_text2.insert(tk.END, "№, Значение в лучшей точке, Лучшая точка \n")
    for i in range(chromo_val):
        line = "{} , {}, {}".format(i,pbest_value[i].round(4), pbest_pos[i].round(4)) + "\n"
        output_text2.insert(tk.END, line)

    resx = "Minimum found: {}".format(gbest_pos.round(6))
    resy = "Minimum value: {}".format(gbest_value.round(6))
    output_text1.insert("1.0", resx + "\n")
    output_text1.insert("1.0", resy + "\n")
    output_text1.insert("1.0", "************************\n")
    print(k)
    if gbest_value <= 10**(-3):
        print("Найденный минимум", gbest_pos)
        print("Значение минимума", gbest_value)
        quit()

def delete():
    output_text1.delete("0.0", tk.END)
    ax.clear()
    canvas.draw()

root = tk.Tk()
root.title("Роевой алгоритм")
root.geometry("1200x800")

left_frame = ttk.Frame(root, padding=10)
right_frame = ttk.LabelFrame(root, padding=10, text="Особи")

root.columnconfigure(0, weight=1)
root.columnconfigure(1, weight=1)

left_frame.grid(row=0, column=0, sticky="nsew")
right_frame.grid(row=0, column=1, sticky="nsew")

frame1 = ttk.LabelFrame(left_frame, padding=10, text="Параметры")
frame11 = ttk.LabelFrame(left_frame, padding=10, text="Количество шагов")
frame2 = ttk.LabelFrame(left_frame, padding=10, text="Управление")
frame3 = ttk.LabelFrame(left_frame, padding=10, text="Результаты")
frame1.pack(side="top", fill="x", padx=5, pady=5)
frame11.pack(side="top", fill="x", padx=5, pady=5)
frame2.pack(side="top", fill="x", padx=5, pady=5)
frame3.pack(side="top", fill="both", expand=True, padx=5, pady=5)

combo_label = ttk.Label(frame1, text="Функция")
combo = ttk.Combobox(frame1, width="40")
combo['value'] = ("8*x[0]**2+4*x[0]*x[1]+5*x[1]**2",)
combo.current(0)

mutant_label = ttk.Label(frame1, text="Размерность")
mutant_entry = ttk.Spinbox(frame1, from_=0, to=10)
mutant_entry.insert(0, "2")
chromoval_label = ttk.Label(frame1, text="Количество особей")
chromoval_entry = ttk.Spinbox(frame1, from_=0, to=10000)
chromoval_entry.insert(0, "50")
minn_label = ttk.Label(frame1, text="Левая граница области поиска")
minn_entry = ttk.Spinbox(frame1, from_=-10000, to=10000)
minn_entry.insert(0, "-50")
maxx_label = ttk.Label(frame1, text="Правая граница области поиска")
maxx_entry = ttk.Spinbox(frame1, from_=-10000, to=10000)
maxx_entry.insert(0, "50")
inercia_label = ttk.Label(frame1, text="Итерационный вес")
inercia_entry = ttk.Spinbox(frame1, from_=-100, to=100)
inercia_entry.insert(0, "0.7")
speed_label = ttk.Label(frame1, text="Коэффициент ускорения")
speed_entry = ttk.Spinbox(frame1, from_=-100, to=100)
speed_entry.insert(0, "1.4")

combo_label.grid(row=0, column=0)
combo.grid(row=1, column=0)

mutant_label.grid(row=2, column=0, padx=0, pady=5)
mutant_entry.grid(row=2, column=1, padx=5, pady=5)

chromoval_label.grid(row=3, column=0, padx=0, pady=5)
chromoval_entry.grid(row=3, column=1, padx=5, pady=5)

minn_label.grid(row=4, column=0, padx=0, pady=5)
minn_entry.grid(row=4, column=1, padx=5, pady=5)

maxx_label.grid(row=5, column=0, padx=0, pady=5)
maxx_entry.grid(row=5, column=1, padx=5, pady=5)

inercia_label.grid(row=6, column=0, padx=0, pady=5)
inercia_entry.grid(row=6, column=1, padx=5, pady=5)

speed_label.grid(row=7, column=0, padx=0, pady=5)
speed_entry.grid(row=7, column=1, padx=5, pady=5)

generations_var = tk.IntVar()
spin = ttk.Spinbox(frame11, from_=10, to=1000, width=10, textvariable=generations_var)
radio_button0 = ttk.Radiobutton(frame11, text="1", variable=generations_var, value=1)
radio_button1 = ttk.Radiobutton(frame11, text="10", variable=generations_var, value=10)
radio_button2 = ttk.Radiobutton(frame11, text="100", variable=generations_var, value=100)
radio_button3 = ttk.Radiobutton(frame11, text="500", variable=generations_var, value=500)
radio_button0.grid(row=0, column=0, padx=1, pady=5)
radio_button1.grid(row=0, column=1, padx=1, pady=5)
radio_button2.grid(row=0, column=2, padx=1, pady=5)
radio_button3.grid(row=0, column=3, padx=1, pady=5)
spin.grid(row=0, column=4, padx=5, pady=5)

k = 0
create_button = ttk.Button(frame2, text="Создать особей", command=create)
calculate_button = ttk.Button(frame2, text="Посчитать шаги", command=main)
k_label = ttk.Label(frame2, text="Пройдено шагов:")
k_entry = ttk.Entry(frame2)
create_button.grid(row=2, column=0, padx=5, pady=5)
calculate_button.grid(row=2, column=1, padx=5, pady=5)
k_label.grid(row=3, column=0, padx=5, pady=5)
k_entry.grid(row=3, column=1, padx=5, pady=5)
output_text1 = ScrolledText(frame3, width=20, height=8)
output_text1.pack(side="top", fill="both", expand=True)
delbutton = ttk.Button(frame3, text="Очистить", command=delete)
delbutton.pack(side="bottom", fill="x")
output_text2 = ScrolledText(right_frame, width=50, height=20)
output_text2.pack()

fig = Figure(figsize=(5, 4), dpi=100)
ax = fig.add_subplot(111)
canvas = FigureCanvasTkAgg(fig, master=right_frame)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

dim = 2
bounds = np.array([int(minn_entry.get()) * dim, int(maxx_entry.get()) * dim])
swarm = np.empty((int(chromoval_entry.get()), dim))
velocity = np.empty((int(chromoval_entry.get()), dim))
pbest_pos = swarm.copy()
pbest_value = np.array([np.inf] * int(chromoval_entry.get()))
gbest_pos = np.zeros(dim)
gbest_value = np.inf
generations_var.set(1)
mystat = []

root.mainloop()