"""
Prog:   dkr1.py
Auth:   Oleksii Krutko, IO-z91
Desc:   TEMK-2. dkr 1. Var 10. 2021
"""
import cmath
import math

from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np

from multiplier_formater import multiple_formatter


def create_axs(n_harmonics):
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
    fig_page1.suptitle("Periodic non-sine currents in linear electrics circuits.",
                       fontsize=16, style='normal')
    fig_page1.subplots_adjust(left=left_margin, right=right_margin, top=top_margin, bottom=bottom_margin,
                              wspace=0.3, hspace=0.2)
    fig_page1.canvas.set_window_title('TEMK-2. Dkr1.')
    axs_page1[0].set_title("Number of harmonics: {}".format(n_harmonics), fontsize=12, color='gray')
    axs_page1[0].set_ylabel(r'$E(\omega t)$', rotation=0, loc='top', fontsize=10, color='gray')
    axs_page1[0].set_xlabel(r'$\omega t$', fontsize=10, loc='right', color='gray')
    axs_page1[0].xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
    axs_page1[0].xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
    axs_page1[0].xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
    axs_page1[0].grid()

    axs_page1[1].set_title('ERC', fontsize=12, color='gray')
    axs_page1[1].set_ylabel(r'$U(\omega t)$', rotation=0, loc='top', fontsize=10, color='gray')
    axs_page1[1].set_xlabel(r'$\omega t$', fontsize=10, style='italic', loc='right', color='gray')
    axs_page1[1].xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
    axs_page1[1].xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
    axs_page1[1].xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
    axs_page1[1].grid()

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
    print("Dkr 1. Var. 10")

    Um = 80
    Rh = 20
    L = 4e-3
    C = 14e-6
    w = 1000.0

    # test
    # Um = 70
    # Rh = 20
    # L = 6e-3
    # C = 14e-6
    # w = 1000.0

    n = 3

    # 1-2
    x_range = generate_x_range(0.0, 2 * math.pi, 1000)
    y_range = generate_y_range(x_range, lambda wt: e(wt, Um))

    axs_page1, fig_page1 = create_axs(n)

    axs_page1[0].plot(x_range, y_range, color='blue')

    Bi = []
    Ci = []

    for i in range(n):
        def b_integrand(wt, Um, harm):
            return ((2 * Um * wt) / math.pi) * math.sin(harm * wt)

        def c_integrand(wt, Um, harm):
            return ((2 * Um * wt) / math.pi) * math.cos(harm * wt)

        Bi.append(quad(b_integrand, 0.0, math.pi / 2, args=(Um, i * 2 + 1))[0] * 2 / math.pi)
        Ci.append(quad(c_integrand, 0.0, math.pi / 2, args=(Um, i * 2 + 1))[0] * 2 / math.pi)

    U = []

    for i in range(n):
        U.append(complex(Bi[i] / math.sqrt(2), Ci[i] / math.sqrt(2)))

    U_m = [math.sqrt(2) * abs(u) for u in U]

    Psi = [cmath.phase(u) for u in U]

    ei_ranges = []

    for i in range(n):
        def ei(wt):
            return U_m[i] * math.sin((2 * i + 1) * wt + Psi[i])

        ei_ranges.append(generate_y_range(x_range, lambda wt: ei(wt)))

        axs_page1[0].plot(x_range, ei_ranges[i], alpha=0.3, linestyle='--', color='red')

    e_sum_range = []

    for i in range(len(x_range)):
        y_sum = 0.0
        for j in range(n):
            y_sum = y_sum + ei_ranges[j][i]

        e_sum_range.append(y_sum)

    axs_page1[0].plot(x_range, e_sum_range, color='green')

    # 3-6

    Um2_ranges = []

    for i in range(n):
        k = i * 2 + 1

        Xl = k * w * L
        Xc = 1 / (k * w * C)

        Um2 = math.sqrt(2) * U[i] * complex(0.0, Xl * Rh) / complex(-Xc ** 2 + 2 * Xc * Xl, Xl * Rh - Xc * Rh)

        Um2_mod = abs(Um2)
        Um2_arg = cmath.phase(Um2)

        def U_m2_i(wt):
            return Um2_mod * math.sin((2 * i + 1) * wt + Um2_arg)

        Um2_ranges.append(generate_y_range(x_range, lambda wt: U_m2_i(wt)))

        axs_page1[1].plot(x_range, Um2_ranges[i], alpha=0.3, linestyle='--', color='red')

    Um2_sum_range = []

    for i in range(len(x_range)):
        y_sum = 0.0
        for j in range(n):
            y_sum = y_sum + Um2_ranges[j][i]

        Um2_sum_range.append(y_sum)

    axs_page1[1].plot(x_range, Um2_sum_range, color='green')

    # 7-8
    print("7-8:")

    k = 1

    Xl = k * w * L
    Xc = 1 / (k * w * C)

    Z1 = complex(0.0, -Xc)
    Z2 = complex(0.0, Xl)
    Z3 = complex(0.0, -Xc)

    A11 = 1 + Z1 / Z2
    A12 = Z1 + Z3 + (Z3 * Z1) / Z2
    A21 = 1 / Z2
    A22 = 1 + Z3 / Z2

    delta = A11 * A22 - A12 * A21

    if abs(delta.real - 1.0) > 1e-6 or abs(delta.imag) > 1e-6:
        print("A-form coefficients wrong!")
        return

    G11 = A21 / A11
    G12 = -delta / A11
    G21 = delta / A11
    G22 = A12 / A11

    print("G11 = {}".format(G11))
    print("G12 = {}".format(G12))
    print("G21 = {}".format(G21))
    print("G22 = {}".format(G22))

    # 9
    print("9:")

    gamma = cmath.log(cmath.sqrt(A11 * A22) + cmath.sqrt(A12 * A21))
    alpha = gamma.real
    beta = gamma.imag
    Zc1 = cmath.sqrt(A12 / A21)

    print("Alpha = {}".format(alpha))
    print("Beta = {}".format(beta))
    print("Zc1 = {}".format(Zc1))

    # 10
    print("10:")

    U2 = 100.0
    print("U2 = {}".format(U2))
    I2 = U2 / Rh
    print("I2 = {}".format(I2))
    U1 = A11 * U2 + A12 * I2
    print("U1 = {}".format(U1))
    I1 = A21 * U2 + A22 * I2
    print("I1 = {}".format(I1))

    plt.show()


main()
