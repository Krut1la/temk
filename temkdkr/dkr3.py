"""
Prog:   dkr3.py
Auth:   Oleksii Krutko, IO-z91
Desc:   TEMK-2. dkr 3. Var 10. 2021
"""
import math

import matplotlib.pyplot as plt

from interpolation import InterpolatorNewton


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
    fig_page1.suptitle("Non-linear circuits with direct current.",
                       fontsize=16, style='normal')
    fig_page1.subplots_adjust(left=left_margin, right=right_margin, top=top_margin, bottom=bottom_margin,
                              wspace=0.3, hspace=0.2)
    fig_page1.canvas.set_window_title('TEMK-2. Dkr3.')
    axs_page1.set_title("", fontsize=12, color='gray')
    axs_page1.set_ylabel(r'$I, A$', rotation=0, loc='top', fontsize=10, color='gray')
    axs_page1.set_xlabel(r'$t, c$', fontsize=10, loc='right', color='gray')
    axs_page1.grid()

    return axs_page1, fig_page1


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
    print("Dkr 3. Var. 10")

    E1 = 75
    E2 = 140
    J = 2
    R1 = 55
    R2 = 60
    R3 = 65
    R4 = 70

    U = (0.0, 35.0, 60.0, 82.0, 100.0, 114.0, 125.0, 144.0, 150.0)
    I = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 4.5)

    I22 = ((J * R2 - E1) * (R1 + R4) + R1 * (E1 + E2))/((R1 + R4)*(R1 + R2 + R3) - R1 ** 2)
    I11 = (E1 + E2 + I22 * R1) / (R1 + R4)

    I1 = I11 - I22
    I2 = J - I22
    I3 = -I22
    I4 = I11
    I5 = I11

    Psource = E1 * I1 + E2 * I5 + J * I2 * R2
    Pcons = (I1 ** 2) * R1 + (I2 ** 2) * R2 + (I3 ** 2) * R3 + (I4 ** 2) * R4

    if math.fabs(Psource - Pcons) > 1e-3:
        print("Currents calculated wrong!")

    Uab_1 = I4 * R4 - I3 * R3
    Uab_2 = I4 * R4 + I1 * R1 + I2 * R2 - E1

    if math.fabs(Uab_2 - Uab_1) > 1e-3:
        print("Currents calculated wrong!")

    R23 = (R2 * R3) / (R1 + R2 + R3)
    R12 = (R1 * R2) / (R1 + R2 + R3)
    R13 = (R1 * R3) / (R1 + R2 + R3)

    Rab = R23 + ((R13 + R4) * R12) / (R4 + R12 + R13)

    axs_page1, fig_page1 = create_axs()

    interp = InterpolatorNewton(U, I)

    x_range = generate_x_range(0.0, Uab_1, 1000)
    y_range_I = generate_y_range(x_range, lambda x: interp.interpolate(x))

    axs_page1.plot(U, I, color='red')

    def I_ne(x):
        return (Uab_1 - x)/Rab

    y_range_Ine = generate_y_range(x_range, lambda x: I_ne(x))

    axs_page1.plot(x_range, y_range_Ine, color='blue')

    for i in range(len(x_range)):
        if math.fabs(y_range_I[i] - y_range_Ine[i]) < 1e-2:
            Une = x_range[i]
            Ine = y_range_Ine[i]
            print("Une = {:.4f}".format(Une))
            print("Ine = {:.4f}".format(Ine))
            axs_page1.axvspan(Une, Une, color="green")
            axs_page1.axhspan(Ine, Ine, color="green")
            break

    plt.show()

    fig_page1.savefig("./input_files/Dkr3_img_png")


main()
