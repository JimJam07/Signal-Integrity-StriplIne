import streamlit as st
from Equation import Expression
from gekko import GEKKO
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from utils import case1, case2, case3

gk = GEKKO()
c = 3 * (10 ** 8)
pi = math.pi
st.title("microwave project Visualisation")
st.sidebar.title("Enter Input Values")

# β = st.sidebar.number_input("enter β value:")
# l = st.sidebar.number_input("enter length of transmission line:")
f = st.sidebar.number_input("enter frequency(Hz)(for left line):")
# f2 = st.sidebar.number_input("enter frequency(Hz)(for right line):")
εr = st.sidebar.number_input("enter relative permitivity:")
eq = st.sidebar.text_input("enter the eq in terms of t:")
fn = Expression(eq, ['t'])
if eq != "":
    t = np.linspace(0, 100, 1000)
    y = []
    for i in t:
        y.append(fn(i))
    figure, ax = plt.subplots(figsize=(7, 3))
    plt.plot(t, y)
    st.pyplot(figure)
type = st.sidebar.selectbox("enter the type:", ("Analysis", "Synthesis"))

if type == "Analysis":
    w = st.sidebar.number_input("enter Width:")
    h = st.sidebar.number_input("enter Height:")
    s = st.sidebar.number_input("enter Spacing:")
    try:
        u = w / h
        g = s / h
        # εe calculation
        a = 1 + 0.02 * (math.log(((u ** 4) + (u / 52) ** 2) / ((u ** 4) + 0.432), math.e)) + math.log(
            1 + ((u / 18.1) ** 3), math.e) / 18.7
        b = 0.564 * (math.pow((εr - 0.9) / (εr + 3), 0.053))
        εe = (εr + 1) / 2 + ((εr - 1) / 2) * (math.pow((1 + (10 / u)), -a * b))
        # Z0 calculation
        F1 = 6 + (2 * math.pi - 6) * (math.exp(-30.666 / u ** 0.7528))
        Z01 = 60 * math.log(F1 / u + math.sqrt(1 + (4 / (u * u))), math.e)
        Z0 = Z01 / math.sqrt(εe)
        st.write(f'Z0 - {Z0}')
        # Even mode
        sharan = 0.8645 * math.pow(u, 0.172)
        si = 1 + g / 1.45 + math.pow(g, 2.09) / 3.95
        alpha = 0.5 * math.exp(-1 * g)
        m_g = 0.2175 + math.pow(4.113 + (math.pow(20.36 / g, 6)), -0.251)
        phi_e = sharan / (si * (alpha * math.pow(u, m_g) + (1 - alpha) * math.pow(u, -1 * m_g)))
        mu = g * (math.exp(-g)) + u * (20 + math.pow(g, 2)) / (10 + math.pow(g, 2))
        eta = 376.73
        Z01e = Z0 / (1 - (Z0 * phi_e / eta))
        F_e = math.pow((1 + 10 / mu), -a * b)
        eee = (εr + 1) / 2 + (εr - 1) * (F_e / 2)
        Z0e = Z01e / (math.sqrt(eee))
        st.write(f'Z0e - {Z0e}')

        p = math.exp(-0.745 * (g ** 0.295)) / math.cosh(g ** 0.68)
        q = math.exp(-1.366 - g)
        r = 1 + 0.15 * (1 - math.exp((1 - (εr - 1) ** 2 / 8.2)) / (1 + g ** (-6)))
        fo1 = 1 - math.exp(-0.179 * (g ** 0.15) - (0.328 * math.pow(g, r)) / math.log(math.e + (g / 7) ** 2.8, math.e))
        fo = fo1 * math.exp(p * math.log(u, math.e) + q * math.sin(math.pi * math.log(u, 10)))
        theta = 1.729 + 1.175 * (math.log(1 + (0.627 / (g + 0.327 * (g ** 2.17))), math.e))
        beta = 0.2306 + (math.log((g ** 10) / (1 + ((g / 3.73) ** 10)), math.e)) / 301.8 + math.log(
            (1 + (0.646 * (g ** 1.175))), math.e) / 5.3
        waste1 = math.log(u, math.e)
        waste2 = math.log(g, math.e)
        n = ((1 / 17.7) + math.exp(-6.424 - 0.76 * waste2 - math.pow(g / 0.23, 5))) * math.log(
            (10 + 68.3 * g * g) / (1 + 32.5 * (g ** 3.093)), math.e)
        # Odd Mode
        phi_o = phi_e - (theta / si) * (math.exp(beta * waste1 * (u ** n)))
        F_o = fo * math.pow((1 + 10 / u), -(a * b))
        eeo = (εr + 1) / 2 + (εr - 1) * (F_o / 2)
        Z01o = Z0 / (1 - (Z0 * phi_o / eta))
        Z0o = Z01o / (math.sqrt(eeo))
        # till this
        st.write(f'Z0o - {Z0o}')
        x = gk.Param(value=np.linspace(1, 20, 100))
        vpe = (c / (math.sqrt(eee)))
        vpo = (c / (math.sqrt(eeo)))
        Le = Z0e / vpe
        Ce = 1 / (Z0e * vpe)
        Lo = Z0o / vpo
        Co = 1 / (Z0o * vpo)

        # C matrix
        c11 = (Ce + Co) / 2
        c12 = (Ce - Co) / 2
        c22 = c11
        c21 = c12

        # L matrix
        l11 = (Le + Lo) / 2
        l12 = (Le - Lo) / 2
        l22 = l11
        l21 = l12
        L = [[l11, l12], [l21, l22]]
        C = [[c11, c12], [c21, c22]]
        st.write(f'L={L}')
        st.write(f'C={C}')
        ω = 2 * pi * f
        # st.write(ω)
        A = -(ω * ω) * (l11 * c11 + l12 * c12)
        B = -(ω * ω) * (l12 * c11 + l11 * c12)
        C = -(ω * ω) * (l21 * c22 + l22 * c21)
        D = -(ω * ω) * (l12 * c12 + l11 * c11)
        noComp = 0
        coeff = [1, 0, -(A + D), 0, (A * D - B * C)]
        roots = np.roots(coeff)
        print(roots)
        A, B = 16.5, 11
        roots2 = [1, -1, 2, -2]
        roots1 = [5, 4]
        roots3 = [6, 2]
        for i in roots:
            if np.iscomplex(i):
                noComp += 1
        case = st.sidebar.selectbox("enter case:", ("case-1", "case-2", "case-3"))
        # case 0
        if case == "case-1":
            posRoots = [t for t in roots2 if t > 0]
            m = posRoots[0]
            n = posRoots[1]
            # print(math.sin(ω))
            time = st.number_input("enter time in sec")
            x, y = case1(m, n, A, time, ω)
            figure, ax = plt.subplots(figsize=(7, 3))
            plt.plot(x, y)
            st.pyplot(figure)
            st.latex(
                r"""\frac{-\left(n^{2}-a\right)\left(\sin wt\right)\cdot e^{-mx}}{\left(m^{2}-n^{2}\right)}+\ \frac{\left(m^{2}-a\right)\left(\sin wt\right)\cdot e^{-nx}}{\left(m^{2}-n^{2}\right)}""")
            st.write(pd.DataFrame({
                "parameter": ["m", "n", "a", "ω"],
                "values": [m, n, A, ω]
            }))

        elif case == "case-2":
            time = st.number_input("enter time(s)")
            v1, v2, x = case2(roots1[0], roots1[1], A, B, time, ω)
            figure, ax = plt.subplots(figsize=(7, 3))
            plt.plot(x, v1)
            st.pyplot(figure)
            figure, ax = plt.subplots(figsize=(7, 3))
            plt.plot(x, v2)
            st.pyplot(figure)
            st.latex(
                r"""V1 = e^{-mx}\left(\left(\sin t\right)\cos nx-\left(\left(n^{2}-m^{2}+a\right)\cdot\frac{\sin t}{2mn}\right)\sin nx\right)""")
            st.latex(
                r"""V2 = \frac{\left(\sin t\right)\left(e^{-mx}\cdot\sin\left(nx\right)\right)}{b}\left(\frac{2a\left(m^{2}-n^{2}\right)-\left(m^{4}+n^{4}\right)-a^{2}-4m^{2}n^{2}}{2mn}\right)""")
            st.write(pd.DataFrame({
                "parameter": ["m", "n", "a", "ω"],
                "values": [roots2[0], roots2[1], A, ω]
            }))
        #     dual mode
        elif case == "case-3":
            time = st.number_input("enter time(s)")
            v1, v2, x = case3(roots3[0], roots3[1], A, B, time, ω)
            print(v1)
            figure, ax = plt.subplots(figsize=(7, 3))
            plt.plot(x, v1)
            st.pyplot(figure)
            figure, ax = plt.subplots(figsize=(7, 3))
            plt.plot(x, v2)
            st.pyplot(figure)
            st.latex(
                r"""V1=\left(\sin\left(2\pi ft\right)\right)\cdot\frac{\left(m^{2}+a\right)}{\left(m^{2}-n^{
                2}\right)}\cos\left(nx\right)\ -\ \left(\sin\left(2\pi ft\right)\right)\cdot\frac{\left(n^{
                2}+a\right)}{\left(m^{2}-n^{2}\right)}\cos\left(mx\right)""")
            st.latex(
                r"""V2= \frac{\left(\sin\left(2\pi ft\right)\right)}{b}\cdot\frac{\left(m^{2}-a\right)\left(n^{
                2}+a\right)}{\left(m^{2}-n^{2}\right)}\cos\left(mx\right)\ -\ \frac{\left(\sin\left(2\pi 
                ft\right)\right)}{b}\cdot\frac{\left(n^{2}-a\right)\left(m^{2}+a\right)}{\left(m^{2}-n^{
                2}\right)}\cos\left(nx\right)""")
            st.write(pd.DataFrame({
                "parameter": ["m", "n", "a", "ω"],
                "values": [roots3[0], roots3[1], A, ω]
            }))
    except ZeroDivisionError:  # to avoid starting error
        pass

    # edge symmetric coupling assumption
