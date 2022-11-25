import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import datetime
import os
import pandas as pd
import scipy.optimize as optimize


def fileparse(nombre_archivo):
    """
    precondicion: recibe un archivo de salida ticra tools

    poscondicion: retorna una tupla con los valores de:
        (frecuencia,
        xmin, xmax, ymin, ymax,
        nfilas, ncolumnas,
        matriz_co, matriz_cx)
    """

    with open(nombre_archivo, encoding="utf8") as file:
        rows = csv.reader(file)

        # sustraccion y transformacion del dato de frecuencia a un float que se va a utilizar luego
        for row in rows:
            if "FREQUENCIES [GHz]:" in row:
                frecuencia = next(rows)
                for c in frecuencia:
                    frecuencia = float(c.strip())
                break

        # se saltean lineas con informacion irrelevante
        for x in range(4):
            next(rows)

        # sustraccion y transformacion de los datos xmin, xmax, ymin e ymax que se van a utilizar luego
        equis = next(rows)
        for c in equis:
            c = c.strip()
            c = " ".join(c.split())
            c = c.split(" ")
            #  hay que llevar primero a float porque es notacion cientifica y no se puede transformar directamente en int
            xmin = int(float(c[0]))
            xmax = int(float(c[2]))
            ymin = int(float(c[1]))
            ymax = int(float(c[3]))

        # sustraccion y transformacion de los datos nfilas y ncolumnas que se van a utilizar luego
        n = next(rows)
        for c in n:
            c = c.strip()
            c = " ".join(c.split())
            c = c.split(" ")
            nfilas = int(c[0])
            ncolumnas = filas = int(c[1])

        # creacion de dos vectores (copolar y crospolar)
        vector_co = np.zeros(0)
        vector_cx = np.zeros(0)

        # sustraccion y movimiento de los datos hacia los vectores
        for fila in rows:
            for c in fila:
                c = c.strip()
                c = " ".join(c.split())
                c = c.split(" ")
                ncomplejo_co = complex(float(c[0]), float(c[1]))
                vector_co = np.append(vector_co, ncomplejo_co)
                ncomplejo_cx = complex(float(c[2]), float(c[3]))
                vector_cx = np.append(vector_cx, ncomplejo_cx)

    # transformacion de los vectores en matrices nfilas, ncolumnas
    matriz_co = vector_co.reshape(nfilas, ncolumnas)
    matriz_cx = vector_cx.reshape(nfilas, ncolumnas)

    return (frecuencia,
            xmin, xmax, ymin, ymax,
            nfilas, ncolumnas,
            matriz_co, matriz_cx)


def procesar(nombre_archivo):
    """
    precondicion: recibe un archivo de salida ticra tools

    poscondicion: llama la funcion fileparse,
    retorna una tupla con los valores de:
        (frecuencia,
         xmin, xmax, ymin, ymax,
         x_range, y_range,
         matriz_co, matriz_cx,
         matriz_co_potencia, matriz_co_pot_norm, matriz_co_pot_norm_log,
         matriz_cx_potencia, matriz_cx_pot_norm, matriz_cx_pot_norm_log,
         E_co_cut, E_co_dB_cut,
         E_cx_cut, E_cx_dB_cut)
    """

    (frecuencia,
     xmin, xmax, ymin, ymax,
     nfilas, ncolumnas,
     matriz_co, matriz_cx) = fileparse(nombre_archivo)

    # creacion de vectores de la grilla de evaluación x e y
    delta_x = (xmax - xmin) / (nfilas - 1)
    x_range = np.arange(xmin, xmax + delta_x, delta_x)
    delta_y = (ymax - ymin) / (ncolumnas - 1)
    y_range = np.arange(ymin, ymax + delta_y, delta_y)

    # creacion de las matrices potencia, normales y logaritmicas a partir de las matrices copolares y crospolares
    matriz_co_potencia = abs(matriz_co)
    maximo_co = np.max(matriz_co_potencia)
    matriz_co_pot_norm = matriz_co_potencia/maximo_co
    # matriz de potencia logaritmica medida en dB
    matriz_co_pot_norm_log = np.log10(matriz_co_pot_norm) * 20

    matriz_cx_potencia = abs(matriz_cx)
    matriz_cx_pot_norm = matriz_cx_potencia/maximo_co
    matriz_cx_pot_norm_log = np.log10(matriz_cx_pot_norm) * 20

    # creacion de los cortes del campo principal (E) y normal (H), a partir de las matrices copolares y crospolares
    E_co_cut = matriz_co_pot_norm[int((nfilas-1)/2)]
    E_cx_cut = matriz_cx_pot_norm[int((nfilas-1)/2)]
    H_co_cut = matriz_co_pot_norm[:, int((nfilas-1)/2)]
    H_cx_cut = matriz_cx_pot_norm[:, int((nfilas-1)/2)]

    E_co_dB_cut = matriz_co_pot_norm_log[int((nfilas-1)/2)]
    E_cx_dB_cut = matriz_cx_pot_norm_log[int((nfilas-1)/2)]
    H_co_dB_cut = matriz_co_pot_norm_log[:, int((nfilas-1)/2)]
    H_cx_dB_cut = matriz_cx_pot_norm_log[:, int((nfilas-1)/2)]

    return(frecuencia,
           xmin, xmax, ymin, ymax,
           nfilas, ncolumnas,
           x_range, y_range,
           matriz_co, matriz_cx,
           matriz_co_potencia, matriz_co_pot_norm, matriz_co_pot_norm_log,
           matriz_cx_potencia, matriz_cx_pot_norm, matriz_cx_pot_norm_log,
           E_co_cut, E_co_dB_cut, H_co_cut, H_co_dB_cut,
           E_cx_cut, E_cx_dB_cut, H_cx_cut, H_cx_dB_cut)


def main(nombre_archivo, Te=12):
    """
    precondicion: recibe un archivo de salida ticra tools y una Te

    poscondicion:  y llama la funcion procesr,
    crea una carpeta llamada "Analysis {fecha y hora} ,
    {filename}", realiza una serie de graficos que guarda en un pdf y,
    calcula los supuestos del modelo y guarda todos los datos relevantes y las
    matrices en un csv llamado "data" dentro de dicha carpeta.
    """

    (frecuencia,
     xmin, xmax, ymin, ymax,
     nfilas, ncolumnas,
     x_range, y_range,
     matriz_co, matriz_cx,
     matriz_co_potencia, matriz_co_pot_norm, matriz_co_pot_norm_log,
     matriz_cx_potencia, matriz_cx_pot_norm, matriz_cx_pot_norm_log,
     E_co_cut, E_co_dB_cut, H_co_cut, H_co_dB_cut,
     E_cx_cut, E_cx_dB_cut, H_cx_cut, H_cx_dB_cut) = procesar(nombre_archivo)

    # calculo de los waistx y waisty expected
    waistx = 0.216 * 8 * (300/frecuencia) * (Te**0.5)
    waisty = waistx

    # creacion matriz X,Y
    X, Y = np.meshgrid(x_range, y_range)

    # creacion de la curva Gaussiana
    Gx = np.exp(-((x_range / waistx) ** 2))  # funcion gaussiana en eje x
    Gx_log = np.log10(np.exp(-((x_range / waistx) ** 2))) * 20  # funcion gaussiana 2D
    G = np.exp(-((X/waistx) ** 2 + (Y/waisty) ** 2))  # funcion gaussiana 2D

    # calculo del Beam coupling
    A = sum(sum(matriz_co * G.conjugate()))
    B = sum(sum(matriz_co * matriz_co.conjugate()))
    C = sum(sum(G * G.conjugate()))
    beam_coupling = ((abs(A)) ** 2) / (B * C).real

    # creacion de ajuste gausiano a los datos propios
    data_xy = []
    for i in range(nfilas):
        for j in range(ncolumnas):
            data_xy.append((x_range[i], y_range[j]))

    z = list(matriz_co_pot_norm.reshape(nfilas * ncolumnas))

    def func(X, amp, x_off, y_off, wx_fit, wy_fit):
        gauss_fit = amp * np.exp(-np.power(X[1] - x_off, 2) / (np.power(wx_fit, 2))-np.power(X[0] - y_off, 2) / (np.power(wy_fit, 2)))
        return gauss_fit

    params, pcov = optimize.curve_fit(func, np.transpose(data_xy), z)
    (amp, x_off, y_off, wx_fit, wy_fit) = params

    # matriz con los datos de la función Gaussiana ajustada (gauss_fit)
    gauss_fit = np.transpose(func([X, Y], amp,
                                  x_off, y_off, wx_fit, wy_fit))

    # Calculo de la gaussicity y elipcicity
    A_fit = sum(sum(matriz_co * gauss_fit.conjugate()))
    B = sum(sum(matriz_co * matriz_co.conjugate()))
    C_fit = sum(sum(gauss_fit * gauss_fit.conjugate()))
    gaussicity = ((abs(A_fit)) ** 2) / (B * C_fit).real
    elipcicity = 1 - (wx_fit / wy_fit)

    # plot Pot co-polar normalizada
    co_polar = plt.figure(figsize=(8, 6), dpi=80)
    plt.pcolor(x_range, y_range, matriz_co_pot_norm,
               cmap='jet',
               vmin=0, vmax=1)
    plt.colorbar()
    plt.grid(linewidth=0.3)
    plt.title(f'E (Co-plar) -normalized- Frequency: {frecuencia} [GHz]')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

    # plot Pot cross-polar normalizada
    cx_polar = plt.figure(figsize=(8, 6), dpi=80)
    plt.pcolor(x_range, y_range, matriz_cx_pot_norm,
               cmap='jet',
               vmin=0, vmax=1)
    plt.colorbar()
    plt.grid(linewidth=0.3)
    plt.title(f'E (Cross-plar) -normalized- Frequency: {frecuencia} [GHz]')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

    # plot co-polar logaritmica normalizada
    co_polar_db = plt.figure(figsize=(8, 6), dpi=80)
    plt.pcolor(x_range, y_range, matriz_co_pot_norm_log,
               cmap='jet',
               vmin=-80, vmax=0)
    plt.colorbar()
    CS_coln = plt.contour(x_range, y_range, matriz_co_pot_norm_log,
                          levels=[-60, -50, -40, -30, -25, -20, -15, -8.7, -3],
                          colors='black')
    plt.clabel(CS_coln, inline=1, fontsize=10)
    plt.grid(linewidth=0.3)
    plt.title(f'E (Co-polar) [dB] -normalized- Frequency: {frecuencia} [GHz]')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

    # plot Pot cross-polar logaritmica normalizada
    cx_polar_db = plt.figure(figsize=(8, 6), dpi=80)
    plt.pcolor(x_range, y_range, matriz_cx_pot_norm_log,
               cmap='jet',
               vmin=-80, vmax=0)
    plt.colorbar()
    CS_cxln = plt.contour(x_range, y_range, matriz_cx_pot_norm_log,
                          levels=[-60, -50, -40, -30, -25, -20, -15, -8.7, -3],
                          colors='black')
    plt.clabel(CS_cxln, inline=1, fontsize=10)
    plt.grid(linewidth=0.3)
    plt.title(f'E (Cross-polar) [dB] -normalized- Frequency: {frecuencia} [GHz]')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

    # plot cortes co-polar y cross polar normalizados
    cut_plot = plt.figure(figsize=(8, 6), dpi=80)
    plt.plot(x_range, E_co_cut,
             color="red",
             linewidth=2.5,
             label="E (co-pol)")
    plt.plot(y_range, H_co_cut,
             color=[0.1, 0.9, 0.5],
             linewidth=2.5,
             label="H (co-pol)")
    plt.plot(x_range, E_cx_cut,
             color=[0, 0.75, 0.75],
             linewidth=2.5,
             label="E (cx-pol)")
    plt.plot(y_range, H_cx_cut,
             color=[0.9290, 0.6940, 0.1250],
             linewidth=2.5,
             label="H (cx-pol)")
    plt.plot(x_range, Gx,
             color="k",
             linewidth=2.5,
             linestyle="--",
             label="Gaussian")
    plt.xlim(xmin, xmax)
    plt.ylim(0, 1.02)
    plt.grid(linewidth=0.3)
    plt.legend(loc='upper right')
    plt.title(f'Power cuts -normalized- Frequency: {frecuencia} [GHz]')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

    # plot cortes co-polar y cross polar logaritmicos normalizados
    dB_cut_plot = plt.figure(figsize=(8, 6), dpi=80)
    plt.plot(x_range, E_co_dB_cut,
             color="red",
             linewidth=2.5,
             label="E (co-pol)")
    plt.plot(y_range, H_co_dB_cut,
             color=[0.1, 0.9, 0.5],
             linewidth=2.5,
             label="H (co-pol)")
    plt.plot(x_range, E_cx_dB_cut,
             color=[0, 0.75, 0.75],
             linewidth=2.5,
             label="E (cx-pol)")
    plt.plot(y_range, H_cx_dB_cut,
             color=[0.9290, 0.6940, 0.1250],
             linewidth=2.5,
             label="H (cx-pol)")
    plt.plot(x_range, Gx_log,
             color="k",
             linewidth=2.5,
             linestyle="--",
             label="Gaussian")
    plt.xlim(xmin, xmax)
    plt.ylim(-80, 1)
    plt.grid(linewidth=0.3)
    plt.legend(loc='upper right')
    plt.title(f'Power cuts [dB] -normalized- Frequency: {frecuencia} [GHz]')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')

    # creacion de carpeta donde se van a guardar las exportaciones
    current_str = datetime.datetime.strftime(datetime.datetime.today(),
                                             "%Y-%m-%d , %H-%M-%S")
    directorio = f'./Analysis {current_str} , {nombre_archivo}'
    os.mkdir(directorio)

    # transformacion df/np a df/pd y creacion de csvs
    index = list(y_range)
    header = list(x_range)
    pdmatriz_co_potencia = pd.DataFrame(matriz_co_potencia)
    pdmatriz_co_potencia.to_csv(f'{directorio}/matriz_co_potencia.csv',
                                index=index, header=header)
    pdmatriz_co_pot_norm = pd.DataFrame(matriz_co_pot_norm)
    pdmatriz_co_pot_norm.to_csv(f'{directorio}/matriz_co_potencia_norm.csv',
                                index=index, header=header)
    pdmatriz_co_pot_norm_log = pd.DataFrame(matriz_co_pot_norm_log)
    pdmatriz_co_pot_norm_log.to_csv(f'{directorio}/matriz_co_potencia_norm_log.csv',
                                    index=index, header=header)
    pdmatriz_cx_potencia = pd.DataFrame(matriz_cx_potencia)
    pdmatriz_cx_potencia.to_csv(f'{directorio}/matriz_cx_potencia.csv',
                                index=index, header=header)
    pdmatriz_cx_pot_norm = pd.DataFrame(matriz_cx_pot_norm)
    pdmatriz_cx_pot_norm.to_csv(f'{directorio}/matriz_cx_potencia_norm.csv',
                                index=index, header=header)
    pdmatriz_cx_pot_norm_log = pd.DataFrame(matriz_cx_pot_norm_log)
    pdmatriz_cx_pot_norm_log.to_csv(f'{directorio}/matriz_cx_potencia_norm_log.csv',
                                    index=index, header=header)
    pdmatriz_gauss_fit = pd.DataFrame(gauss_fit)
    pdmatriz_gauss_fit.to_csv(f'{directorio}/matriz_gauss_fit.csv',
                              index=index, header=header)

    # guardado de plots en pdfs
    with PdfPages(f'{directorio}/figures.pdf') as pdf:
        pdf.savefig(co_polar)
        pdf.savefig(cx_polar)
        pdf.savefig(co_polar_db)
        pdf.savefig(cx_polar_db)
        pdf.savefig(cut_plot)
        pdf.savefig(dB_cut_plot)

    # creacion del csv data con info reelevante
    to_write = [
        [f'filename,{nombre_archivo}'],
        [f'frecuencia,{frecuencia}'],
        [f'xmin,{xmin}'],
        [f'xmax,{xmax},'],
        [f'ymin,{ymin}'],
        [f'ymax,{ymax}'],
        [f'nfilas,{nfilas}'],
        [f'ncolumnas,{ncolumnas}'],
        [''],
        ['ajuste gaussiano'],
        [f'amp,{amp}'],
        [f'x_off,{x_off}'],
        [f'y_off,{y_off}'],
        [f'waistx_expected,{waistx}'],
        [f'waisty_expected,{waisty}'],
        [f'waistx_fit,{wx_fit}'],
        [f'waisty_fit,{wy_fit}'],
        [''],
        [f'elipcicity,{elipcicity}'],
        [f'gaussicity,{gaussicity}'],
        [f'beam_coupling,{beam_coupling}'],
        ]
    with open(f'{directorio}/data.csv', encoding="utf8",
              mode="w", newline="") as data:
        writer = csv.writer(data, delimiter=",")
        writer.writerows(to_write)


main("CASS_planar_grid_B2_3_Casseg_MF")


if __name__ == "__main__":
    argvs = sys.argv
    if len(argvs) == 3:
        nombre_archivo = argvs[1]
        Te = argvs[2]
        main(nombre_archivo, Te)
    else:
        print('Argumentos necesarios: NombreScript ruta')
