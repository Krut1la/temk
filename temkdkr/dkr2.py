"""
Prog:   dkr2.py
Auth:   Oleksii Krutko, IO-z91
Desc:   TEMK-2. dkr 2. Var 10. 2021
"""
import cmath
import math

from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np

from multiplier_formater import multiple_formatter


def create_axs():
    # A4 page
    fig_width_cm = 21
    fig_height_cm = 29.7

    extra_margin_cm = 1
    margin_left_cm = 2 + extra_margin_cm
    margin_right_cm = 1.0
    margin_bottom_cm = 1.0 + extra_margin_cm
    margin_top_cm = 3.0 + extra_margin_cm
    inches_per_cm = 1 / 2.54

    fig_width = fig_width_cm * inches_per_cm  # width in inches
    fig_height = fig_height_cm * inches_per_cm  # height in inches

    margin_left = margin_left_cm * inches_per_cm
    margin_right = margin_right_cm * inches_per_cm
    margin_bottom = margin_bottom_cm * inches_per_cm
    margin_top = margin_top_cm * inches_per_cm

    left_margin = margin_left / fig_width
    right_margin = 1 - margin_right / fig_width
    bottom_margin = margin_bottom / fig_height
    top_margin = 1 - margin_top / fig_height

    fig_size = [fig_width, fig_height]

    plt.rc('figure', figsize=fig_size)

    fig_page1, axs_page1 = plt.subplots(1, 1)
    fig_page1.suptitle("Transient response.",
                       fontsize=16, style='normal')
    fig_page1.subplots_adjust(left=left_margin, right=right_margin, top=top_margin, bottom=bottom_margin,
                              wspace=0.3, hspace=0.2)
    fig_page1.canvas.set_window_title('TEMK-2. Dkr2.')
    axs_page1.set_title("", fontsize=12, color='gray')
    axs_page1.set_ylabel(r'$I, A$', rotation=0, loc='top', fontsize=10, color='gray')
    axs_page1.set_xlabel(r'$t, c$', fontsize=10, loc='right', color='gray')
    axs_page1.grid()

    return axs_page1, fig_page1


def e(wt, Um):
    if 0.0 <= wt < math.pi / 2:
        return 2 * Um * wt / math.pi

    if math.pi / 2 <= wt < math.pi:
        return 0.0

    if math.pi <= wt < 3 * math.pi / 2:
        return -2 * Um * (wt - math.pi) / math.pi

    if 3 * math.pi / 2 <= wt < 2 * math.pi:
        return 0.0


def generate_x_range(a, b, m):
    x = []
    for i in range(0, m + 1):
        x.append(a + i * (b - a) / m)

    return x


def generate_y_range(x_range, func):
    y = []
    for x in x_range:
        y.append(func(x))

    return y


def main():
    print("Dkr 2. Var. 10")

    U = 65
    R1 = 100
    R2 = 100
    R3 = 100
    L = 0.19
    C = 50e-6
    T_r = 1.5

    # test
    # U = 20
    # R1 = 40
    # R2 = 20
    # R3 = 40
    # L = 1e-2
    # # C = 50e-6
    # T_r = 4e-3

    axs_page1, fig_page1 = create_axs()

    p = -(R2*R3 + R1*R2 + R1*R3)/(L*(R1 + R3))

    A = 1/(R1 + R3)

    print("p = {:.04}".format(p))
    print("A = {:.04}".format(A))

    def u_t(t):
        if 0 <= t <= T_r:
            return t*U/T_r
        if T_r <= t <= 1000000.0:
            return 0.0

    def gL_t_x(t):
        return A*(1 - math.exp(p*(t)))

    def integrand_1(x, t):
        return (U/T_r)*gL_t_x(t - x)

    def integrand_2(x, t):
        return (U/T_r)*gL_t_x(t - x)

    def iL_t_1(t):
        return u_t(0.0)*gL_t_x(t) + quad(integrand_1, 0.0, t, args=t)[0]

    def iL_t_2(t):
        return u_t(0.0)*gL_t_x(t) + quad(integrand_2, 0.0, T_r, args=t)[0] - u_t(T_r)*gL_t_x(t - T_r)

    x_range = generate_x_range(0.0, T_r, 1000)
    y_range = generate_y_range(x_range, lambda t: iL_t_1(t))

    axs_page1.plot(x_range, y_range, color='green')

    x_range = generate_x_range(T_r, 6e-3, 1000)
    y_range = generate_y_range(x_range, lambda t: iL_t_2(t))

    axs_page1.plot(x_range, y_range, color='green')

    plt.show()

    # fig_page1.savefig("./input_files/Dkr1_img_png")


main()
