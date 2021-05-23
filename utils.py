from Equation import Expression
from gekko import GEKKO
import matplotlib.pyplot as plt
import math
import numpy as np


def case1(m, n, A, t, w):
    x = np.linspace(0, 10, 1000)
    y = ((m * m - A) * math.sin(w * t) * np.exp(-n * x)) / (m * m - n * n) - (
                (n * n - A) * math.sin(w * t) * np.exp(-m * x)) / (m * m - n * n)
    return x, y


def case2(m, n, A, B, t, w):
    x = np.linspace(0, 10, 1000)
    v1 = np.exp(-m * x) * (
                (math.sin(w*t) * np.cos(n * x)) - ((n * n - m * m + A) * math.sin(w*t) / (2 * m * n)) * np.sin(n * x))
    # v2 calculation
    temp = (2 * A * (m * m - n * n) - (m ** 4 + n ** 4) - A * A - 4 * m * m * n * n) / (2 * m * n)
    v2 = (math.sin(w*t) * (np.exp(-m * x) * np.sin(n * x)) * temp) / B
    return v1, v2, x


def case3(m, n, A, B, t, w):
    x = np.linspace(0, 10, 1000)
    v1 = math.sin(w * t) * ((m * m + A) / (m * m - n * n)) * np.cos(n * x) - math.sin(w * t) * (
                (n * n + A) / (m * m - n * n)) * np.cos(m * x)
    temp = math.sin(w * t)/B
    temp1 = ((m * m - A) * (n * n + A) * np.cos(m * x) / (m * m - n * n)) - (
                (n * n - A) * (m * m + A) * np.cos(m * x) / (m * m - n * n))
    v2 = temp1 * temp
    print("in")
    return v1, v2, x
