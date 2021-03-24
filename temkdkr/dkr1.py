"""
Prog:   dkr1.py
Auth:   Oleksii Krutko, IO-z91
Desc:   TEMK-2. dkr 1. Var 10. 2021
"""
import math

from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np


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

    fig_page1, axs_page1 = plt.subplots(2, 1)
    fig_page1.suptitle("test",
                       fontsize=16, style='normal')
    fig_page1.subplots_adjust(left=left_margin, right=right_margin, top=top_margin, bottom=bottom_margin,
                              wspace=0.3, hspace=0.2)
    fig_page1.canvas.set_window_title('Lab 3. Polynomial interpolation. Var. 10')
    axs_page1[0].set_title("Test", fontsize=12, color='gray')
    axs_page1[0].set_ylabel(r'$f(x)$', rotation=0, loc='top', fontsize=10, color='gray')
    axs_page1[0].set_xlabel(r'$x$', fontsize=10, loc='right', color='gray')
    axs_page1[0].grid()

    axs_page1[1].set_title('Nodes, interpolated', fontsize=12, color='gray')
    axs_page1[1].set_ylabel("f(x)", rotation=0, loc='top', fontsize=10, color='gray')
    axs_page1[1].set_xlabel(r'$x$', fontsize=10, style='italic', loc='right', color='gray')
    axs_page1[1].grid()

    return axs_page1, fig_page1

def e(wt, Um):
    if 0.0 <= wt < math.pi/2:
        return 2 * Um * wt / math.pi

    if math.pi/2 <= wt < math.pi:
        return 0.0

    if math.pi <= wt < 3*math.pi/2:
        return -2 * Um * (wt - math.pi) / math.pi

    if 3*math.pi/2 <= wt < 2*math.pi:
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
    print("Dkr 1. Var. 10")

    # Um = 80
    # Rh = 20
    # L = 4e-3
    # C = 14e-6
    # w = 1000.0

    # test
    Um = 70
    Rh = 20
    L = 6e-3
    C = 14e-6
    w = 1000.0

    # a = 0.0
    #
    # b = math.pi/2
    #
    # I = quad(integrand, a, b, args=Um)[0]*2/math.pi

    # draw original func
    x_range = generate_x_range(0.0, 2*math.pi, 1000)
    y_range = generate_y_range(x_range, lambda wt: e(wt, Um))

    axs_page1, fig_page1 = create_axs()

    axs_page1[0].plot(x_range, y_range, color='blue')


    n = 3

    B = []
    C = []

    for i in range(n):
        def b_integrand(wt, Um, harm):
            return ((2 * Um * wt) / math.pi) * math.sin(harm * wt)

        def c_integrand(wt, Um, harm):
            return ((2 * Um * wt) / math.pi) * math.cos(harm * wt)

        B.append(quad(b_integrand, 0.0, math.pi / 2, args=(Um, i * 2 + 1))[0] * 2 / math.pi)
        C.append(quad(c_integrand, 0.0, math.pi / 2, args=(Um, i * 2 + 1))[0] * 2 / math.pi)

    U = []

    for i in range(n):
        U.append(complex(B[i]/math.sqrt(2), C[i]/math.sqrt(2)))

    Um = [math.sqrt(2)*abs(u) for u in U]


    print(U)
    # plt.show()


main()